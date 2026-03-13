Security Policy
===============

Scope
-----

This project is designed for local workstation usage by default.

- Default runtime is localhost-only.
- Remote access is disabled unless explicitly enabled via `INITBINDER_ALLOW_REMOTE=true`.

Threat Model (Local Mode)
-------------------------

Expected:

- Single-user local machine
- Trusted local filesystem
- User-managed cluster credentials and SSH setup outside the app

Not guaranteed in local mode:

- Multi-user internet-facing deployment hardening
- Strong authentication/authorization for public network exposure

If You Enable Remote Access
---------------------------

If you set `INITBINDER_ALLOW_REMOTE=true`, you are responsible for:

- Network isolation / firewall policy
- Reverse proxy and TLS termination
- Authentication and access control
- Resource controls (rate limits, quotas)

Reporting a Vulnerability
-------------------------

Please report suspected vulnerabilities privately to the maintainers instead of opening a public issue with exploit details.

When reporting, include:

- Affected endpoint or workflow
- Reproduction steps
- Impact assessment
- Suggested mitigation (if known)
