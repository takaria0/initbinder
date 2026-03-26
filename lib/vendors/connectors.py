from __future__ import annotations

import asyncio
import csv
import dataclasses
import json
import re
import time
from dataclasses import dataclass, asdict
from typing import Any, Dict, Iterable, List, Optional

import httpx
from bs4 import BeautifulSoup
from fake_useragent import UserAgent
from tenacity import retry, wait_exponential, stop_after_attempt

# ------------------------------
# Common types & utilities
# ------------------------------

@dataclass
class AntigenRecord:
    vendor: str
    gene_symbol: Optional[str]
    protein_name: Optional[str]
    sku: str
    url: str
    species: List[str]
    expression_host: Optional[str]
    sequence: Optional[str]
    tags: List[str]
    has_biotin: bool
    has_his: bool
    has_fc: bool
    has_avi: bool
    macs_ready: bool

    def to_row(self) -> Dict[str, Any]:
        d = asdict(self)
        # Flatten list fields for CSV friendliness
        d["species"] = ",".join(self.species)
        d["tags"] = ",".join(self.tags)
        return d


def _norm_whitespace(t: str) -> str:
    return re.sub(r"\s+", " ", (t or "").strip())


# --- replace this function entirely ---
def _map_tags(raw_tags: Iterable[str]) -> List[str]:
    """Normalize vendor tag strings to a small vocabulary without conflating Avi with Biotin."""
    out: set[str] = set()
    for t in raw_tags:
        key = t.lower().replace(" ", "").replace("__", "_")
        parts = re.split(r"[-_]+", key)
        for p in parts:
            if p in {"his", "h", "his6", "his6x"}:
                out.add("His")
            elif p in {"fc", "hfc", "rfc", "mfc", "llama", "igg1", "igg4"}:
                out.add("Fc")
            elif p in {"avi", "avi_tag"}:
                out.add("Avi")
            elif p in {"biotin", "biotinylated"}:
                out.add("Biotin")
            elif p in {"flag"}:
                out.add("FLAG")
            elif p in {"gst"}:
                out.add("GST")
            elif p in {"mbp"}:
                out.add("MBP")
            elif p in {"strep", "strep_ii"}:
                out.add("Strep")
            elif p in {"myc"}:
                out.add("Myc")
    # NOTE: Do NOT infer Biotin from Avi anymore.
    return sorted(out)

# --- replace this function entirely ---
def _flags_from_tags(tags: List[str]) -> Dict[str, bool]:
    s = set(tags)
    return {
        "has_biotin": "Biotin" in s,
        "has_his": "His" in s,
        "has_fc": "Fc" in s,
        "has_avi": "Avi" in s,
        # macs_ready is now STRICT: requires explicit Biotin
        "macs_ready": "Biotin" in s,
    }


# ------------------------------
# Base connector
# ------------------------------

class BaseVendorConnector:
    def __init__(self, *, mode: str = "headless", timeout: float = 25.0, throttle_s: float = 1.2):
        # accepted backends:
        #   - 'interactive' => Playwright visible window
        #   - 'headless'    => Playwright headless (no window; good for scale)
        #   - 'xhr'         => static HTTP fetch + HTML parse (fallback to headless if JS required)
        assert mode in {"interactive", "headless", "xhr"}
        self.mode = mode
        self.timeout = timeout
        self.ua = UserAgent()
        self.throttle_s = throttle_s

    def _sleep(self):
        time.sleep(self.throttle_s)

    # Public API
    def search_proteins(
        self,
        query: str,
        *,
        species_preference: Optional[Iterable[str]] = None,
        limit: int = 60,
    ) -> List[Dict[str, Any]]:
        raise NotImplementedError

    # CSV helper
    @staticmethod
    def write_csv(rows: List[Dict[str, Any]], path: str):
        if not rows:
            return
        keys = list(rows[0].keys())
        with open(path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=keys)
            writer.writeheader()
            for r in rows:
                writer.writerow(r)


# ------------------------------
# Sino Biological connector
# ------------------------------
class SinoBioConnector(BaseVendorConnector):
    BASE = "https://www.sinobiological.com"
    # SEARCH_URL_TPL = BASE + "/search?keywords={kw}&searchCode=1"  # 1 = Protein
    # https://www.sinobiological.com/search/by-category?keywords=PD1&categoryCode=1&searchType=1
    SEARCH_URL_TPL = BASE + "/search/by-category?keywords={kw}&categoryCode=1&searchType=1"

    def _browser_impl(
        self, query: str, species_preference: Optional[Iterable[str]], limit: int, *, headless: bool
    ) -> List[Dict[str, Any]]:
        try:
            from playwright.sync_api import sync_playwright, TimeoutError as PWTimeout
        except ImportError as e:
            raise RuntimeError(
                "Playwright is required for browser mode. Install with:\n"
                "  pip install playwright && python -m playwright install chromium"
            ) from e

        results: List[AntigenRecord] = []
        species_pref = {s.lower() for s in (species_preference or [])}
        debug_dump = bool(getattr(self, "debug_dump", False))
        # slow motion only when interactive; speed up headless
        slow_mo = int(getattr(self, "slow_mo", 120 if not headless else 0))

        from urllib.parse import quote

        html = ""
        with sync_playwright() as pw:
            browser = pw.chromium.launch(headless=headless, slow_mo=slow_mo)
            ua = getattr(self, "ua", None)
            ua_str = getattr(ua, "chrome", None) or (
                "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
                "AppleWebKit/537.36 (KHTML, like Gecko) "
                "Chrome/124.0.0.0 Safari/537.36"
            )
            context = browser.new_context(
                user_agent=ua_str,
                locale="en-US",
                timezone_id="America/Los_Angeles",
                viewport={"width": 1366, "height": 900},
                java_script_enabled=True,
            )
            context.add_init_script("Object.defineProperty(navigator, 'webdriver', {get: () => undefined});")
            page = context.new_page()

            try:
                search_url = self.SEARCH_URL_TPL.format(kw=quote(query))
                page.goto(search_url, wait_until="domcontentloaded", timeout=int(self.timeout * 1000))
                try:
                    page.wait_for_load_state("networkidle", timeout=8000)
                except Exception:
                    pass

                page.wait_for_selector(
                    "#search_result_cntr, #search_result_body",
                    timeout=int(self.timeout * 1000),
                )
                for sel in ("#mask_loading", ".loading", ".spinner", ".lds-ellipsis"):
                    try:
                        page.wait_for_selector(sel, state="hidden", timeout=4000)
                    except Exception:
                        pass

                page.wait_for_function(
                    """
                    () => {
                    const hasRows = document.querySelectorAll('#search_result_cntr ul.result-item.body').length > 0
                                || document.querySelectorAll('#search_result_cntr .result-item.body').length > 0;
                    const bodyTxt = (document.querySelector('#search_result_body')?.innerText || '');
                    const noRes = /no results|not\\s*found|0\\s*results/i.test(bodyTxt);
                    return hasRows || noRes;
                    }
                    """,
                    timeout=int(self.timeout * 1000),
                )

                for _ in range(10):
                    rows_now = len(
                        page.query_selector_all("#search_result_cntr ul.result-item.body, #search_result_cntr .result-item.body")
                    )
                    if rows_now >= limit:
                        break
                    page.mouse.wheel(0, 1600)
                    self._sleep()

                html = page.content()

            except PWTimeout:
                html = page.content()
                if debug_dump:
                    try:
                        with open("sino_debug_timeout.html", "w", encoding="utf-8") as f:
                            f.write(html)
                    except Exception:
                        pass
                raise
            finally:
                # keep a tiny pause for interactive viewing; no-op in headless
                try:
                    if not headless:
                        page.wait_for_timeout(600)
                except Exception:
                    pass
                context.close()
                browser.close()

        # --- parse (unchanged from your original _visible_impl) ---
        soup = BeautifulSoup(html, "lxml")
        cards = soup.select("#search_result_cntr ul.result-item.body[data-cate='protein']")
        if not cards:
            candidates = soup.select("#search_result_cntr ul.result-item.body, #search_result_cntr .result-item.body")
            cards = [ul for ul in candidates if ul.select_one("li[data-col='species']") and ul.select_one("li[data-col='expressionhost']")]

        parsed: List[Dict[str, Any]] = []
        for ul in cards[:limit]:
            try:
                sku_a = ul.select_one("li[data-col='catalogue'] a")
                sku = _norm_whitespace(sku_a.get_text()) if sku_a else ""
                url = (sku_a.get("href") or "") if sku_a else ""
                if url.startswith("/"):
                    url = self.BASE + url

                desc_li = ul.select_one("li[data-col='description']")
                protein_name = _norm_whitespace(desc_li.get_text()) if desc_li else None
                description_full_text = (desc_li.get_text(" ", strip=True) if desc_li else "").lower()

                sp_li = ul.select_one("li[data-col='species']")
                species_txt = _norm_whitespace(sp_li.get_text()) if sp_li else ""
                species = [s.strip() for s in species_txt.split(",") if s.strip()]

                host_li = ul.select_one("li[data-col='expressionhost']")
                expression_host = _norm_whitespace(host_li.get_text()) if host_li else None

                seq_li = ul.select_one("li[data-col='sequence']")
                sequence = _norm_whitespace(seq_li.get_text()) if seq_li else None

                raw_tag_keys = []
                for span in ul.select("li[data-col='tag'] span.lbl"):
                    key = span.get("data-fixpos", "")
                    if "|" in key:
                        parts = key.split("|")
                        if len(parts) >= 2:
                            raw_tag_keys.append(parts[1])
                if "biotin" in description_full_text:
                    raw_tag_keys.append("biotin")

                tags = _map_tags(raw_tag_keys)
                flags = _flags_from_tags(tags)

                gene_symbol = None
                if protein_name:
                    m = re.search(r"\b([A-Z0-9]{2,7})(?:/|\b)", protein_name)
                    if m:
                        gene_symbol = m.group(1)

                if species_pref and species and not any(s.lower() in species_pref for s in species):
                    continue

                rec = AntigenRecord(
                    vendor="Sino Biological",
                    gene_symbol=gene_symbol,
                    protein_name=protein_name,
                    sku=sku,
                    url=url,
                    species=species,
                    expression_host=expression_host,
                    sequence=sequence,
                    tags=tags,
                    **flags,
                )
                parsed.append(rec.to_row())
            except Exception:
                continue

        return parsed

    # --- XHR backend (kept as-is; will auto-upgrade to visible if JS rendering is required)
    @retry(wait=wait_exponential(multiplier=1.2, min=1, max=8), stop=stop_after_attempt(3))
    def _xhr_impl(
        self, query: str, species_preference: Optional[Iterable[str]], limit: int
    ) -> List[Dict[str, Any]]:
        headers = {"User-Agent": getattr(self.ua, "chrome", "Mozilla/5.0"), "Accept-Language": "en-US,en;q=0.9"}
        params = {"categoryCode": "", "keywords": query}
        with httpx.Client(headers=headers, timeout=self.timeout, follow_redirects=True) as client:
            r = client.get(self.BASE + "/search/by-category", params=params)
            r.raise_for_status()
            html = r.text

        soup = BeautifulSoup(html, "lxml")
        # if not soup.select("#search_result_cntr ul.result-item.body"):
        #     # requires JS; switch to visible browser
        #     return self._visible_impl(query, species_preference, limit)
        # return self._parse_sino_rendered(soup, limit, species_preference)
        if not soup.select("#search_result_cntr ul.result-item.body"):
            # requires JS; switch to headless browser (no window)
            return self._browser_impl(query, species_preference, limit, headless=True)
        return self._parse_sino_rendered(soup, limit, species_preference)


    @staticmethod
    def _parse_sino_rendered(
        soup: BeautifulSoup, limit: int, species_preference: Optional[Iterable[str]]
    ) -> List[Dict[str, Any]]:
        species_pref = {s.lower() for s in (species_preference or [])}
        out: List[Dict[str, Any]] = []
        cards = soup.select("#search_result_cntr ul.result-item.body[data-cate='protein']")
        for ul in cards[:limit]:
            sku_a = ul.select_one("li[data-col='catalogue'] a")
            sku = _norm_whitespace(sku_a.get_text()) if sku_a else ""
            url = sku_a["href"] if sku_a and sku_a.has_attr("href") else ""
            if url and url.startswith("/"):
                url = SinoBioConnector.BASE + url

            desc_li = ul.select_one("li[data-col='description']")
            protein_name = _norm_whitespace(desc_li.get_text()) if desc_li else None
            # --- MODIFIED: Check full description text for biotin ---
            description_full_text = (desc_li.get_text(" ", strip=True) if desc_li else "").lower()

            sp_li = ul.select_one("li[data-col='species']")
            species_txt = _norm_whitespace(sp_li.get_text()) if sp_li else ""
            species = [s.strip() for s in species_txt.split(",") if s.strip()]
            if species_pref and species and not any(s.lower() in species_pref for s in species):
                continue

            host_li = ul.select_one("li[data-col='expressionhost']")
            expression_host = _norm_whitespace(host_li.get_text()) if host_li else None

            seq_li = ul.select_one("li[data-col='sequence']")
            sequence = _norm_whitespace(seq_li.get_text()) if seq_li else None

            raw_tag_keys = []
            for span in ul.select("li[data-col='tag'] span.lbl"):
                key = span.get("data-fixpos", "")
                if "|" in key:
                    parts = key.split("|")
                    if len(parts) >= 2:
                        raw_tag_keys.append(parts[1])

            # --- MODIFIED: Add biotin to tags if found in description ---
            if "biotin" in description_full_text:
                raw_tag_keys.append("biotin")
            
            tags = _map_tags(raw_tag_keys)
            flags = _flags_from_tags(tags)

            gene_symbol = None
            if protein_name:
                m = re.search(r"\b([A-Z0-9]{2,7})(?:/|\b)", protein_name)
                if m:
                    gene_symbol = m.group(1)

            rec = AntigenRecord(
                vendor="Sino Biological",
                gene_symbol=gene_symbol,
                protein_name=protein_name,
                sku=sku,
                url=url,
                species=species,
                expression_host=expression_host,
                sequence=sequence,
                tags=tags,
                **flags,
            )
            out.append(rec.to_row())
        return out

    # public
    # def search_proteins(
    #     self,
    #     query: str,
    #     *,
    #     species_preference: Optional[Iterable[str]] = None,
    #     limit: int = 60,
    # ) -> List[Dict[str, Any]]:
    #     """
    #     For this vendor we default to **visible browser** so you can inspect behavior.
    #     Even if self.mode == 'headless', we will still use the visible browser.
    #     Use self.mode == 'xhr' for the lightweight parser (it auto-upgrades to visible if needed).
    #     """
    #     backend = (getattr(self, "mode", "") or "visible").lower()
    #     if backend == "xhr":
    #         return self._xhr_impl(query, species_preference, limit)
    #     # Treat anything else (including 'headless') as visible browser:
    #     return self._visible_impl(query, species_preference, limit)


    def search_proteins(
        self,
        query: str,
        *,
        species_preference: Optional[Iterable[str]] = None,
        limit: int = 60,
    ) -> List[Dict[str, Any]]:
        """
        Backends:
        - mode='xhr'         : static HTTP; auto-fallback to headless if JS needed
        - mode='headless'    : Playwright headless (no visible window)
        - mode='interactive' : Playwright with a visible window (debug)
        """
        backend = (getattr(self, "mode", "") or "headless").lower()
        if backend == "xhr":
            return self._xhr_impl(query, species_preference, limit)
        elif backend == "headless":
            return self._browser_impl(query, species_preference, limit, headless=True)
        else:
            return self._browser_impl(query, species_preference, limit, headless=False)


# ------------------------------
# ACROBiosystems Connector (NEW)
# ------------------------------
import re
import httpx
from typing import Any, Dict, Iterable, List, Optional
from bs4 import BeautifulSoup
from urllib.parse import quote

# Assumptions:
# - BaseVendorConnector provides: self.ua.chrome, self.timeout, self.mode, self._sleep(), self.throttle_s
# - AntigenRecord dataclass
# - Helpers: _norm_whitespace, _map_tags, _flags_from_tags
#
# Install Playwright first if needed:
#   pip install playwright && python -m playwright install chromium

class ACROConnector(BaseVendorConnector):
    BASES = [
        "https://jp.acrobiosystems.com",
        "https://www.acrobiosystems.com",
    ]

    # Current search URL shapes observed:
    #   /search-keywords-<slug>.html?is_search=1
    #   /search.php?keywords=<query>&is_search=1
    # (older legacy)
    #   /search?key=<query>
    SEARCH_PATTERNS = [
        "/search-keywords-{slug}.html?is_search=1",
        "/search.php?keywords={q}&is_search=1",
        "/search?key={q}",
    ]

    # ---------- Visible browser (Playwright) ----------
    def _visible_impl(
        self, query: str, species_preference: Optional[Iterable[str]], limit: int
    ) -> List[Dict[str, Any]]:
        try:
            from playwright.sync_api import sync_playwright, TimeoutError as PWTimeout
        except ImportError as e:
            raise RuntimeError(
                "Playwright is required for browser mode. Install with:\n"
                "  pip install playwright && python -m playwright install chromium"
            ) from e

        species_pref = {s.lower() for s in (species_preference or [])}
        debug_dump = bool(getattr(self, "debug_dump", False))
        slow_mo = int(getattr(self, "slow_mo", 120))  # slower = easier to see

        html = ""
        with sync_playwright() as pw:
            browser = pw.chromium.launch(headless=False, slow_mo=slow_mo)

            ua_str = getattr(getattr(self, "ua", None), "chrome", None) or (
                "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
                "AppleWebKit/537.36 (KHTML, like Gecko) "
                "Chrome/124.0.0.0 Safari/537.36"
            )
            context = browser.new_context(
                user_agent=ua_str,
                locale="en-US",
                timezone_id="America/Los_Angeles",
                viewport={"width": 1366, "height": 900},
                java_script_enabled=True,
            )
            # Hide webdriver flag for friendlier behavior
            context.add_init_script("Object.defineProperty(navigator, 'webdriver', {get: () => undefined});")
            page = context.new_page()

            # Build candidate URLs (JP + WWW; slug + query variants)
            slug = re.sub(r"\s+", "-", query.strip())
            q = quote(query.strip())
            candidates = []
            for base in self.BASES:
                for pat in self.SEARCH_PATTERNS:
                    u = base + pat.format(slug=quote(slug), q=q)
                    if u not in candidates:
                        candidates.append(u)

            try:
                usable = False
                for u in candidates:
                    page.goto(u, wait_until="domcontentloaded", timeout=int(self.timeout * 1000))
                    # networkidle can hang; best-effort
                    try:
                        page.wait_for_load_state("networkidle", timeout=8000)
                    except Exception:
                        pass

                    # Wait until we either see a results table or a "no results" message in the right column
                    try:
                        page.wait_for_selector("div.right table, .search_result_right table", timeout=int(self.timeout * 1000))
                        usable = True
                        break
                    except Exception:
                        # try next candidate
                        continue

                if not usable:
                    # Capture page for debugging
                    html = page.content()
                    if debug_dump:
                        try:
                            with open("acro_debug_no_table.html", "w", encoding="utf-8") as f:
                                f.write(html)
                        except Exception:
                            pass
                    # Give up
                    return []

                # Scroll to encourage any lazy rendering; pause between scrolls
                for _ in range(10):
                    # Stop early if we already see enough rows
                    if len(page.query_selector_all("div.right table tbody tr, .search_result_right table tbody tr")) >= limit:
                        break
                    page.mouse.wheel(0, 1600)
                    self._sleep()

                # Visually outline the table for human inspection
                try:
                    page.evaluate("""
                        const t = document.querySelector('div.right table, .search_result_right table');
                        if (t) { t.style.outline = '3px solid #e91e63'; t.style.scrollMargin = '40px'; t.scrollIntoView({behavior: 'smooth', block: 'center'}); }
                    """)
                except Exception:
                    pass

                # Snapshot final HTML for parsing
                html = page.content()

            except PWTimeout:
                html = page.content()
                if debug_dump:
                    try:
                        with open("acro_debug_timeout.html", "w", encoding="utf-8") as f:
                            f.write(html)
                    except Exception:
                        pass
                raise
            finally:
                # leave the window up a moment before closing
                try:
                    page.wait_for_timeout(600)
                except Exception:
                    pass
                context.close()
                browser.close()

        # Parse rendered HTML
        soup = BeautifulSoup(html, "lxml")
        return self._parse_acro_rendered(soup, limit, species_preference)

    # ---------- XHR (static fetch) ----------
    @retry(wait=wait_exponential(multiplier=1.2, min=1, max=8), stop=stop_after_attempt(3))
    def _xhr_impl(
        self, query: str, species_preference: Optional[Iterable[str]], limit: int
    ) -> List[Dict[str, Any]]:
        headers = {"User-Agent": getattr(self.ua, "chrome", "Mozilla/5.0"), "Accept-Language": "en-US,en;q=0.9"}
        slug = re.sub(r"\s+", "-", query.strip())
        q = quote(query.strip())

        # try JP first (more consistent markup), then WWW
        for base in self.BASES:
            for pat in self.SEARCH_PATTERNS:
                url = base + pat.format(slug=quote(slug), q=q)
                try:
                    with httpx.Client(headers=headers, timeout=self.timeout, follow_redirects=True) as client:
                        r = client.get(url)
                        r.raise_for_status()
                        soup = BeautifulSoup(r.text, "lxml")
                        rows = self._parse_acro_rendered(soup, limit, species_preference)
                        if rows:
                            return rows
                except Exception:
                    continue

        # If static didn’t yield a parseable table, open a visible browser
        return self._visible_impl(query, species_preference, limit)

    # ---------- Parser for rendered HTML ----------
    @staticmethod
    def _build_header_index(tbl) -> Dict[str, int]:
        # Accept JP + EN
        synonyms = {
            "molecule": {"分子", "molecule"},
            "catalog": {"製造番号", "catalog", "cat", "cat no", "cat. no."},
            "species": {"種類", "species"},
            "description": {"製品説明", "description"},
            "host": {"ホスト", "host"},
            "structure": {"構造", "structure"},
            "purity": {"純度", "purity"},
            "feature": {"特徴", "feature", "features"},
            "price": {"規格と価格", "price"},
            # "compare": {"比べる", "compare"}  # not used
        }

        # ACRO’s header uses <td> in <thead>
        cells = tbl.select("thead tr:first-child > *")
        heads = [ _norm_whitespace(c.get_text()) for c in cells ]
        heads_l = [ h.lower() for h in heads ]

        idx: Dict[str, int] = {}
        for key, names in synonyms.items():
            for name in names:
                if name in heads_l:
                    idx[key] = heads_l.index(name)
                    break
        return idx

    @staticmethod
    def _extract_tags_from_text(text: str) -> List[str]:
        """
        Extract tag-ish tokens from description/feature cells.
        Examples: Biotinylated, His, (Mouse IgG2a) Fc, Avitag, MALS verified, HPLC-verified, PE.
        """
        t = (text or "").lower()
        tokens = []
        if "biotin" in t:
            tokens.append("Biotinylated")
        if re.search(r"\bfc\b|igg", t):
            tokens.append("Fc")
        if re.search(r"\bhis\b|his tag", t):
            tokens.append("His")
        if "avi" in t:  # Avitag / AviTag
            tokens.append("Avitag")
        if "mals" in t:
            tokens.append("MALS verified")
        if "hplc" in t:
            tokens.append("HPLC verified")
        if re.search(r"\bpe\b|pe-", t):
            tokens.append("PE")
        return _map_tags(tokens)

    def _parse_acro_rendered(
        self, soup: BeautifulSoup, limit: int, species_preference: Optional[Iterable[str]]
    ) -> List[Dict[str, Any]]:
        species_pref = {s.lower() for s in (species_preference or [])}

        # Prefer the main right-column results table
        tbl = soup.select_one("div.right table") or soup.select_one(".search_result_right table") or soup.select_one("table")
        if not tbl:
            return []

        header_idx = self._build_header_index(tbl)
        if "catalog" not in header_idx or "description" not in header_idx:
            return []

        out: List[Dict[str, Any]] = []
        for tr in tbl.select("tbody tr"):
            tds = tr.find_all("td")
            if len(tds) < max(header_idx.values()) + 1:
                continue
            try:
                # SKU + product URL (usually both present as anchors)
                cat_cell = tds[header_idx["catalog"]]
                desc_cell = tds[header_idx["description"]]

                sku_a = cat_cell.select_one("a[href]") or desc_cell.select_one("a[href]")
                sku = _norm_whitespace(sku_a.get_text()) if sku_a else _norm_whitespace(cat_cell.get_text())
                if not sku:
                    continue

                href = ""
                if sku_a and sku_a.has_attr("href"):
                    href = sku_a["href"]
                else:
                    # very rare: pull from the description link
                    a_alt = desc_cell.select_one("a[href]")
                    if a_alt and a_alt.has_attr("href"):
                        href = a_alt["href"]

                if not href:
                    # last resort: look for compare button with URL in onclick
                    cmp_a = tr.select_one("a.compare-btn[onclick*=\"/P\"]")
                    if cmp_a:
                        m = re.search(r"'(\/P[^']+\.html)'", cmp_a.get("onclick", ""))
                        if m:
                            href = m.group(1)

                if not href:
                    continue

                if href.startswith("/"):
                    # use canonical domain if present; JP base otherwise
                    base = soup.select_one("link[rel='canonical']")
                    if base and base.has_attr("href"):
                        m = re.match(r"https?://[^/]+", base["href"])
                        base_domain = m.group(0) if m else self.BASES[0]
                    else:
                        base_domain = self.BASES[0]
                    url = base_domain + href
                else:
                    url = href

                # Species
                species = []
                if "species" in header_idx:
                    sp_txt = _norm_whitespace(tds[header_idx["species"]].get_text(" "))
                    if sp_txt:
                        species = [s.strip() for s in re.split(r"[,/;、]", sp_txt) if s.strip()]

                # Host
                expression_host = None
                if "host" in header_idx:
                    host_txt = _norm_whitespace(tds[header_idx["host"]].get_text(" "))
                    expression_host = host_txt or None

                # Molecule / protein name (e.g., "PD-1")
                protein_name = None
                if "molecule" in header_idx:
                    protein_name = _norm_whitespace(tds[header_idx["molecule"]].get_text()) or None

                # Tags from description + (optional) structure/purity/feature cells
                desc_txt = _norm_whitespace(desc_cell.get_text(" "))
                extra_txt = []
                for k in ("structure", "purity", "feature"):
                    if k in header_idx:
                        extra_txt.append(_norm_whitespace(tds[header_idx[k]].get_text(" ")))
                tags = self._extract_tags_from_text(" | ".join([desc_txt] + extra_txt))
                flags = _flags_from_tags(tags)

                # Gene symbol heuristic: "... / PDCD1 Protein ..." or molecule cell
                gene_symbol = None
                m = re.search(r"/\s*([A-Z0-9\-]{2,15})\s+Protein", desc_txt)
                if m:
                    gene_symbol = m.group(1)
                elif protein_name and re.fullmatch(r"[A-Z0-9\-]{2,15}", protein_name):
                    gene_symbol = protein_name

                # species filter
                if species_pref and species and not any(s.lower() in species_pref for s in species):
                    continue

                out.append(
                    AntigenRecord(
                        vendor="ACROBiosystems",
                        gene_symbol=gene_symbol,
                        protein_name=protein_name,
                        sku=sku,
                        url=url,
                        species=species,
                        expression_host=expression_host,
                        sequence=None,
                        tags=tags,
                        **flags,
                    ).to_row()
                )
                if len(out) >= limit:
                    break
            except Exception:
                continue

        return out

    # ---------- Public API ----------
    def search_proteins(
        self,
        query: str,
        *,
        species_preference: Optional[Iterable[str]] = None,
        limit: int = 60,
    ) -> List[Dict[str, Any]]:
        """
        For ACRO we also default to a **visible browser** so you can inspect behavior.
        Set mode='xhr' for a lightweight fetch (it falls back to visible if necessary).
        """
        backend = (getattr(self, "mode", "") or "visible").lower()
        if backend == "xhr":
            return self._xhr_impl(query, species_preference, limit)
        # Treat anything else (including 'headless') as visible:
        return self._visible_impl(query, species_preference, limit)


# ------------------------------
# Small CLI for quick checks
# ------------------------------

def _print_preview(rows: List[Dict[str, Any]], n: int = 5):
    for r in rows[:n]:
        print(
            f"{r['vendor']:<18} {r['sku']:<16} tags={r['tags']:<25} species={r['species']:<15} url={r['url']}"
        )


def main():
    import argparse
    ap = argparse.ArgumentParser("Harvest purchasable antigen options from vendors")
    ap.add_argument("query", help="Protein/gene keyword, e.g., 'PD-1' or 'EGFR'.")
    ap.add_argument("--vendor", choices=["sino", "acro", "both"], default="acro")
    # ap.add_argument("--mode", choices=["headless", "xhr"], default="headless")
    ap.add_argument("--mode", choices=["interactive", "headless", "xhr"], default="headless")
    ap.add_argument("--species", nargs="*", default=None)
    ap.add_argument("--limit", type=int, default=60)
    ap.add_argument("--out", help="Optional path to CSV output")
    args = ap.parse_args()

    all_rows: List[Dict[str, Any]] = []
    if args.vendor in {"sino", "both"}:
        sino = SinoBioConnector(mode=args.mode)
        rows = sino.search_proteins(args.query, species_preference=args.species, limit=args.limit)
        print(f"Sino Biological: {len(rows)} candidates")
        _print_preview(rows)
        all_rows.extend(rows)
    if args.vendor in {"acro", "both"}:
        acro = ACROConnector(mode=args.mode)
        rows = acro.search_proteins(args.query, species_preference=args.species, limit=args.limit)
        print(f"ACROBiosystems: {len(rows)} candidates")
        _print_preview(rows)
        all_rows.extend(rows)

    if args.out and all_rows:
        BaseVendorConnector.write_csv(all_rows, args.out)
        print(f"Wrote {len(all_rows)} rows -> {args.out}")


if __name__ == "__main__":
    main()
