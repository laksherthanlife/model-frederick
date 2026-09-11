# Williamson 2009 -- Complete cAMP Model of *S. cerevisiae*

Williamson T, Schwartz J-M, Kell DB, Stateva L (2009). "Deterministic mathematical models
of the cAMP pathway in *Saccharomyces cerevisiae*." *BMC Systems Biology* **3**:70.
doi:10.1186/1752-0509-3-70 -- PMID 19607691 -- PMC2719611.

## Licence

Creative Commons Attribution 2.0 Generic (CC BY 2.0),
<https://creativecommons.org/licenses/by/2.0/>, as declared in the article's own
`<permissions>` block: "Copyright (c) 2009 Williamson et al; licensee BioMed Central Ltd.
This is an Open Access article distributed under the terms of the Creative Commons
Attribution License, which permits unrestricted use, distribution, and reproduction in any
medium, provided the original work is properly cited."

Redistribution here is therefore permitted, with the attribution above. This is a stronger
position than `jalihal2021/`, which carries no licence and is custody-only and gitignored.

## Files

| file | sha256 | what it is |
|---|---|---|
| `1752-0509-3-70-S1.xml` | `e608d6651dd3332c0a2f146bdd9402c77ee07f7c6358284bef911da5bcd1b782` | Additional file 1: the Complete cAMP Model, SBML Level 2 Version 1 |
| `PMC2719611.xml` | `73769eb3c3cb45b415de3ede58dbad3082191cf90788313cb3cc7ad28a13e58d` | article full text (NCBI efetch, `db=pmc&id=2719611&retmode=xml`), carrying Tables 1-5 and the Figure 8/9 captions |

The SBML was retrieved from the Europe PMC supplementary-file archive
(`https://www.ebi.ac.uk/europepmc/webservices/rest/PMC2719611/supplementaryFiles`); the
publisher's own `static-content.springer.com` copy answers HTTP 403 to an anonymous client.

## Why this directory exists

`src/ystwin/mech/carbon_williamson2009.py` is a hand transcription of this SBML. The shared
loader `mech/kinetic_sbml.py` **refuses** the file at line 332 ("rate and algebraic rules are
unsupported"), because it carries seven `rateRule`s. Those rules are audited rather than
stripped -- see `RATE_RULE_AUDIT` in that module.
