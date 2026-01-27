# Security Policy

## Supported Versions

This project is provided primarily as a research/academic tool for running HingeProt in Google Colab/Jupyter environments.  
Security fixes are applied to the **latest `main` branch**. Older tags/commits may not receive security updates.

| Version/Branch | Supported          |
| -------------- | ------------------ |
| `main`         | :white_check_mark: |
| Older tags     | :x:                |

> If you are using this project in a production-like environment, we strongly recommend pinning dependencies and regularly updating to the latest commit on `main`.

---

## Reporting a Vulnerability

If you believe you have found a security vulnerability, please **do not open a public GitHub issue**.

Instead, report it privately using one of the following methods:

1. **GitHub Security Advisories (preferred)**  
   Go to the repository’s **Security** tab → **Report a vulnerability**.

2. **Email**  
   Send details to: **enesemretas@gmail.com**  

### What to include
Please include as much of the following as possible:
- A clear description of the issue and potential impact
- Steps to reproduce (proof-of-concept is helpful)
- Affected files/components (e.g., `ui.py`, Colab upload logic, shell/Perl execution)
- Environment details (Colab vs local Jupyter, OS, Python version)
- Any relevant logs or screenshots (avoid sharing secrets)

### Response timeline
- **Acknowledgement:** within **2 days**
- **Status update:** within **7 days** after acknowledgement, if more time is needed

### Disclosure
If the report is accepted, we will work on a fix and coordinate a responsible disclosure timeline.  
If the report is declined (e.g., not a security issue or out of scope), we will explain the reasoning.

---

## Scope Notes (Important)

Because this project runs external tools (Perl scripts, compiled binaries) and may execute shell commands in notebook environments:

- Treat notebook environments as **trusted** (especially when running in Colab with your own account).
- Avoid running untrusted code or untrusted PDB files from unknown sources.
- Do not store secrets (API keys/tokens) in notebooks or in the repository.

Thank you for helping keep the project safe.
