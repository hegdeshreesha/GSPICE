# Security policy

GSPICE 1.3 is an academic beta. Only the newest tagged beta receives security
fixes. Netlists and model libraries are untrusted input; run unknown decks and
OSDI binaries in an isolated environment.

Please report vulnerabilities privately to the repository owner through
GitHub's private vulnerability reporting feature. Include the affected version,
platform, minimal reproducer, impact, and any suggested mitigation. Do not
attach confidential PDK or foundry data.

OSDI libraries are native executable code loaded into the simulator process.
GSPICE does not sandbox them or verify their signatures.
