#!/usr/bin/env python3
"""
Parse ACROBiosystems search results HTML and emit a TSV compatible with
targets_catalog/webscraper/acrobio_biotinylated_unique.tsv.

Usage examples:
  python targets_catalog/webscraper/acrobio_biotin_pipeline.py --html "/path/to/downloads/biotinylated _ ACROBiosystems.html"
  python targets_catalog/webscraper/acrobio_biotin_pipeline.py --url "https://www.acrobiosystems.com/search?keywords=biotinylated"
  python targets_catalog/webscraper/acrobio_biotin_pipeline.py --url "https://www.acrobiosystems.com/search?keywords=biotinylated" --page_param page
  python targets_catalog/webscraper/acrobio_biotin_pipeline.py --url "https://www.acrobiosystems.com/search?keywords=biotinylated" --mode playwright --headed
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from urllib.parse import parse_qs, urlencode, urlparse, urlunparse, urljoin

import pandas as pd
import requests
from bs4 import BeautifulSoup

OUTPUT_COLUMNS = [
    "Catalog",
    "Description",
    "Species",
    "Expression Host",
    "Sequence",
    "Tag",
    "Activity",
    "Images",
    "sku",
    "url",
    "protein_name",
    "gene_symbol",
    "tags",
    "species",
    "expression_host",
    "target_name",
    "antigen_url",
]


def _read_html(url: Optional[str], html_path: Optional[str]) -> str:
    if html_path:
        return Path(html_path).read_text(errors="ignore")
    if not url:
        raise ValueError("Provide --url or --html")
    headers = {
        "User-Agent": (
            "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
            "AppleWebKit/537.36 (KHTML, like Gecko) "
            "Chrome/120.0.0.0 Safari/537.36"
        )
    }
    resp = requests.get(url, headers=headers, timeout=30)
    resp.raise_for_status()
    return resp.text


def _read_html_playwright(url: str, wait_selector: str, timeout_ms: int, headed: bool) -> str:
    try:
        from playwright.sync_api import sync_playwright  # type: ignore
    except Exception as exc:  # pragma: no cover - optional dependency
        raise RuntimeError("Playwright is required for --mode playwright.") from exc

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=not headed)
        page = browser.new_page()
        page.goto(url, wait_until="domcontentloaded")
        if wait_selector:
            page.wait_for_selector(wait_selector, timeout=timeout_ms)
        html = page.content()
        browser.close()
    return html


def scrape_acrobio_playwright_paginated(
    url: str,
    max_pages: int,
    wait_selector: str,
    timeout_ms: int,
    headed: bool,
) -> Tuple[List[Dict[str, str]], int]:
    try:
        from playwright.sync_api import sync_playwright  # type: ignore
    except Exception as exc:  # pragma: no cover - optional dependency
        raise RuntimeError("Playwright is required for --mode playwright.") from exc

    all_rows: List[Dict[str, str]] = []
    seen: set[str] = set()
    page_count = 0

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=not headed)
        page = browser.new_page()
        page.goto(url, wait_until="domcontentloaded")
        if wait_selector:
            page.wait_for_selector(wait_selector, timeout=timeout_ms)

        for _ in range(max_pages):
            page_count += 1
            html = page.content()
            try:
                page_rows = parse_acrobio_table(html)
            except RuntimeError:
                break

            new_count = 0
            for row in page_rows:
                key = _unique_key(row)
                if key not in seen:
                    seen.add(key)
                    all_rows.append(row)
                    new_count += 1

            next_btn = page.query_selector("button.btn-next")
            if not next_btn:
                break
            aria_disabled = (next_btn.get_attribute("aria-disabled") or "").lower()
            if aria_disabled == "true":
                break

            try:
                prev_first = page.eval_on_selector(
                    "tr.el-table__row td",
                    "el => el.innerText.trim()",
                )
            except Exception:
                prev_first = None

            next_btn.click()

            if prev_first:
                try:
                    page.wait_for_function(
                        "prev => { const el = document.querySelector('tr.el-table__row td'); return el && el.innerText.trim() !== prev; }",
                        prev_first,
                        timeout=timeout_ms,
                    )
                except Exception:
                    page.wait_for_timeout(1500)
            else:
                page.wait_for_timeout(1500)

            if new_count == 0 and page_count > 1:
                break

        browser.close()

    return all_rows, page_count


def _normalize_header(text: str) -> str:
    return "".join(ch.lower() for ch in text if ch.isalnum())


def _find_results_table(soup: BeautifulSoup):
    required = {
        "molecule",
        "catno",
        "species",
        "host",
        "productdescription",
    }
    for table in soup.find_all("table"):
        headers = [th.get_text(" ", strip=True) for th in table.find_all("th")]
        header_norm = {_normalize_header(h) for h in headers}
        if required.issubset(header_norm):
            return table
    # fallback: try any tbody with el-table__row
    tbody = soup.find("tbody")
    if tbody and tbody.find("tr", class_="el-table__row"):
        return tbody
    return None


def _cell_text(cell) -> str:
    return cell.get_text(" ", strip=True) if cell else ""


def _extract_tags(cat_cell, sku: str) -> str:
    label = cat_cell.find(class_="label") if cat_cell else None
    if label:
        tag_text = " ".join(label.stripped_strings)
        return tag_text.strip()
    # fallback: all text in cell except sku
    parts = [p.strip() for p in cat_cell.stripped_strings] if cat_cell else []
    parts = [p for p in parts if p and p != sku]
    return " ".join(parts).strip()


def _extract_sku_and_url(cat_cell, base_url: str) -> (str, str):
    if not cat_cell:
        return "", ""
    link = cat_cell.find("a", href=True)
    if link:
        href = link["href"]
        return link.get_text(" ", strip=True), urljoin(base_url, href)
    return cat_cell.get_text(" ", strip=True), ""


def parse_acrobio_table(html: str, base_url: str = "https://www.acrobiosystems.com") -> List[Dict[str, str]]:
    soup = BeautifulSoup(html, "lxml")
    table = _find_results_table(soup)
    rows = []
    if table:
        rows = table.find_all("tr", class_="el-table__row")
        if not rows:
            rows = soup.find_all("tr", class_="el-table__row")
        if not rows:
            rows = table.find_all("tr")
    if not rows:
        rows = soup.find_all("tr", class_="el-table__row")
    if not rows:
        raise RuntimeError("Could not locate the ACROBiosystems results table rows in the HTML.")

    parsed: List[Dict[str, str]] = []
    for row in rows:
        tds = row.find_all("td", recursive=False)
        if len(tds) < 5:
            continue

        target_name = _cell_text(tds[0])
        sku, url = _extract_sku_and_url(tds[1], base_url)
        tags = _extract_tags(tds[1], sku)
        species = _cell_text(tds[2])
        expression_host = _cell_text(tds[3])
        description = _cell_text(tds[4])
        structure = _cell_text(tds[5]) if len(tds) > 5 else ""
        purity = _cell_text(tds[6]) if len(tds) > 6 else ""
        feature = _cell_text(tds[7]) if len(tds) > 7 else ""

        # Consolidate tag-like fields
        tag_field = " ".join([t for t in [tags, feature, purity] if t]).strip()
        if "biotin" in description.lower() and "biotin" not in tag_field.lower():
            tag_field = (tag_field + " Biotinylated").strip()

        antigen_url = url or (f"https://www.acrobiosystems.com/search?keywords={sku}" if sku else "")

        parsed.append(
            {
                "Catalog": sku,
                "Description": description,
                "Species": species,
                "Expression Host": expression_host,
                "Sequence": structure,
                "Tag": tag_field,
                "Activity": feature,
                "Images": "",
                "sku": sku,
                "url": url,
                "protein_name": description,
                "gene_symbol": target_name,
                "tags": tag_field,
                "species": species,
                "expression_host": expression_host,
                "target_name": target_name,
                "antigen_url": antigen_url,
            }
        )

    return parsed


def _extract_max_page(html: str) -> Optional[int]:
    soup = BeautifulSoup(html, "lxml")
    numbers: List[int] = []
    for li in soup.select(".el-pager li.number"):
        text = li.get_text(" ", strip=True)
        if text.isdigit():
            numbers.append(int(text))
    if numbers:
        return max(numbers)
    return None


def _build_page_url(
    base_url: str,
    page: int,
    page_param: str,
    page_size: Optional[int] = None,
    page_size_param: str = "pageSize",
) -> str:
    parsed = urlparse(base_url)
    query = parse_qs(parsed.query)
    query[page_param] = [str(page)]
    if page_size is not None:
        query[page_size_param] = [str(page_size)]
    new_query = urlencode(query, doseq=True)
    return urlunparse(parsed._replace(query=new_query))


def _unique_key(row: Dict[str, str]) -> str:
    for key in ("sku", "url", "Catalog", "Description"):
        val = (row.get(key) or "").strip()
        if val:
            return val
    return str(row)


def scrape_acrobio_paginated(
    url: str,
    max_pages: int,
    page_param: str,
    page_size: Optional[int],
    page_size_param: str,
    no_paginate: bool,
    mode: str,
    wait_selector: str,
    timeout_ms: int,
    headed: bool,
) -> Tuple[List[Dict[str, str]], int]:
    all_rows: List[Dict[str, str]] = []
    seen: set[str] = set()
    page = 1
    empty_pages = 0
    total_pages: Optional[int] = None

    while True:
        if mode == "playwright":
            if no_paginate:
                html = _read_html_playwright(url, wait_selector, timeout_ms, headed)
                page_rows = parse_acrobio_table(html)
                all_rows = page_rows
                return all_rows, 1
            return scrape_acrobio_playwright_paginated(
                url,
                max_pages=max_pages,
                wait_selector=wait_selector,
                timeout_ms=timeout_ms,
                headed=headed,
            )

        page_url = _build_page_url(url, page, page_param, page_size, page_size_param) if not no_paginate else url
        html = _read_html(page_url, None)
        page_rows = parse_acrobio_table(html)

        if total_pages is None:
            total_pages = _extract_max_page(html)

        new_count = 0
        for row in page_rows:
            key = _unique_key(row)
            if key not in seen:
                seen.add(key)
                all_rows.append(row)
                new_count += 1

        if not page_rows or new_count == 0:
            empty_pages += 1
        else:
            empty_pages = 0

        if no_paginate:
            break
        if total_pages and page >= total_pages:
            break
        if page >= max_pages:
            break
        if empty_pages >= 2:
            break

        page += 1

    return all_rows, page


def build_dataframe(rows: List[Dict[str, str]], dedupe: bool) -> pd.DataFrame:
    df = pd.DataFrame(rows)
    for col in OUTPUT_COLUMNS:
        if col not in df.columns:
            df[col] = ""
    df = df[OUTPUT_COLUMNS]

    if dedupe:
        df = df[df["target_name"].astype(str).str.strip() != ""]
        df = df.drop_duplicates(subset=["target_name"]).sort_values("target_name")
    return df


def main() -> None:
    default_root = Path(__file__).resolve().parent
    default_out = default_root / "acrobio_biotinylated_unique.tsv"
    default_raw = default_root / "acrobio_biotinylated_raw.csv"
    parser = argparse.ArgumentParser(description="Parse ACROBiosystems biotinylated search HTML.")
    parser.add_argument("--url", type=str, default=None, help="Search URL to fetch.")
    parser.add_argument("--html", type=str, default=None, help="Path to saved HTML file.")
    parser.add_argument("--output", type=str, default=str(default_out), help="TSV output path.")
    parser.add_argument("--raw_output", type=str, default=str(default_raw), help="Raw CSV output path.")
    parser.add_argument("--no_dedupe", action="store_true", help="Disable dedupe by target_name.")
    parser.add_argument("--max_pages", type=int, default=200, help="Maximum pages to fetch when using --url.")
    parser.add_argument("--page_param", type=str, default="page", help="Query parameter name for page index.")
    parser.add_argument("--page_size", type=int, default=None, help="Page size value to request (optional).")
    parser.add_argument("--page_size_param", type=str, default="pageSize", help="Query parameter name for page size.")
    parser.add_argument("--no_paginate", action="store_true", help="Disable pagination even when --url is provided.")
    parser.add_argument(
        "--mode",
        type=str,
        choices=["requests", "playwright"],
        default="requests",
        help="Fetch mode for URLs (use playwright for JS-rendered pages).",
    )
    parser.add_argument(
        "--wait_selector",
        type=str,
        default=".el-table__row",
        help="CSS selector to wait for in playwright mode.",
    )
    parser.add_argument(
        "--timeout_ms",
        type=int,
        default=30000,
        help="Playwright wait timeout in milliseconds.",
    )
    parser.add_argument(
        "--headed",
        action="store_true",
        help="Use a visible browser window in playwright mode.",
    )
    args = parser.parse_args()

    if args.html and not args.url:
        html = _read_html(None, args.html)
        rows = parse_acrobio_table(html)
        if not rows:
            raise RuntimeError("No rows parsed from ACROBiosystems HTML.")
        fetched_pages = 1
    else:
        if not args.url:
            raise ValueError("Provide --url or --html")
        rows, fetched_pages = scrape_acrobio_paginated(
            args.url,
            max_pages=args.max_pages,
            page_param=args.page_param,
            page_size=args.page_size,
            page_size_param=args.page_size_param,
            no_paginate=args.no_paginate,
            mode=args.mode,
            wait_selector=args.wait_selector,
            timeout_ms=args.timeout_ms,
            headed=args.headed,
        )
        if not rows:
            raise RuntimeError("No rows parsed from ACROBiosystems HTML.")

    df_raw = pd.DataFrame(rows)
    df = build_dataframe(rows, dedupe=not args.no_dedupe)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)
    raw_path = Path(args.raw_output)
    raw_path.parent.mkdir(parents=True, exist_ok=True)
    df_raw.to_csv(raw_path, index=False)

    print(f"[ok] parsed {len(rows)} rows across {fetched_pages} page(s) -> {len(df)} unique targets")
    print(f"[ok] wrote {out_path}")
    print(f"[ok] wrote {raw_path}")


if __name__ == "__main__":
    main()
