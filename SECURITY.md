# Security Policy

## Supported versions

| Version | Supported |
|---|---|
| 1.0.x | Yes |
| < 1.0 | No |

## Reporting a vulnerability

Please report security issues privately through
[GitHub Security Advisories](https://github.com/bekirhankarakurt/radar-signal-processing-sim/security/advisories/new)
instead of opening a public issue. Include the version or commit, the scenario file or
command line that triggers the problem and what you observed.

You can expect an acknowledgement within 7 days and a fix or mitigation plan within 30 days
for confirmed issues.

## Scope

`radar_sim` reads scenario files and writes output files to a directory you choose. Relevant
reports include memory-safety bugs triggered by malformed scenario files, path handling
problems and denial of service through crafted input. The simulator processes synthetic data
only and contains no networking code.
