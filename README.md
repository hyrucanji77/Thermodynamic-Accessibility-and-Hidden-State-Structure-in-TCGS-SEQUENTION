# Thermodynamic Accessibility and Hidden-State Structure in TCGS–SEQUENTION

**Henry Arellano-Peña · manuscript version 4 · 2 November 2026**

This repository accompanies *Thermodynamic Accessibility and Hidden-State Structure in TCGS–SEQUENTION: Environmental Memory, Epistemic Information, and an Exact Finite-Time Work Law*.

## Availability

The executable numerical package is complete: the three original Python files, four original source workbooks, thirteen result files (including both verification reports), and four generated figures are present. The source workbooks and the two publisher supplementary PDFs were accepted only after exact SHA-256 verification.

**Eleven PDFs are not yet mirrored:** the compiled manuscript (`main.pdf`), nine individual framework references, and the exact experimental-article PDF used by the manuscript. These files are supplied in the companion PDF distribution and identified in `COMPANION_PDFS.sha256`. The publisher download did not match the recorded article checksum and was not silently substituted. A relative PDF link in the reference index may therefore refer to a companion-distribution file rather than an existing repository file.

LaTeX sources, Overleaf archives, compilation containers, earlier manuscript drafts, and private editorial records are excluded.

## Scientific scope

The paper derives the exact finite-time accessible work for the specified passive two-coordinate harmonic model with a free final trap centre. It distinguishes the total conditional free-energy bound, the hidden conditional resource, and work accessible under stated observation and control restrictions. The mathematical results are conditional on the declared model; the TCGS–SEQUENTION interpretation is distinguished from the conventional operational explanation.

The unchanged version 4 manuscript has 40 pages, 27 numbered references, four figures, and seven tables. The numerical implementation also evaluates the retrospective, no-branch-fit comparison with the publisher workbooks of Muruga and colleagues. This is not independent experimental reconstruction from raw trajectories.

## Reproduce the results

Python 3.13.5 and the pinned dependencies in `requirements.txt` reproduce the calculations:

```sh
python -m pip install -r requirements.txt
python code/reproduce.py
python code/label_resolution.py
```

The original workbooks are included; the calculations do not require network access. The optional command `python tools/fetch_publisher_inputs.py` retrieves missing publisher inputs only when their bytes match the recorded checksums.

| File | Role |
| --- | --- |
| `code/accessibility.py` | Riccati solution, optimal mean path, constant-speed work, Gaussian conditioning, and independent discrete-control optimization. |
| `code/reproduce.py` | Original calculations, three figures, workbook summaries, and 70 checks. |
| `code/label_resolution.py` | Ramp comparisons, all 24 label assignments, conditional-mean magnitude inversion, score decomposition, initial-jump prediction, fourth figure, and 35 checks. |

**70/70 and 35/35 checks pass, totaling 105.** Both suites were also executed successfully in [GitHub Actions](https://github.com/hyrucanji77/Thermodynamic-Accessibility-and-Hidden-State-Structure-in-TCGS-SEQUENTION/actions/runs/41544265228). These are mathematical and numerical consistency checks, not 105 independent experimental tests. The code uses nanometres, seconds, and zeptojoules. No branch-specific parameter fitting is performed.

## Workbook interpretation

The selected assignment is **B → 11, D → 00, F → 01, H → 10**; the original workbook labels and values are unchanged. The denominator score identifies D, H, and the B/F pair; full work curves distinguish B from F. Both scores reuse the same model and data and are not independent statistical confirmations.

Under transferred stiffness and temperature, the D baseline residual corresponds to a +4.80320 nm conditional-mean magnitude offset, compared with +5.40195 nm for F. Baseline inversion determines magnitude, not sign. Finite bandwidth is a possible, not established, explanation for the residuals. This comparison is retrospective parameter transfer, not parameter-free physics, prospective validation, or source-author confirmation of the metadata assignment.

## References and attribution

[REFERENCE_INDEX.html](REFERENCE_INDEX.html) and [references.json](references.json) preserve all 27 numbered references and their version notes. `REFERENCE_MANIFEST.sha256` identifies the exact twelve reference PDFs and four original workbooks in the complete companion distribution. Framework papers are identified individually.

The exact 71-page *No Final Map* copy defines SR in §2.4, p.5. Its equivalence to PhilArchive version 1, deposited 15 August 2026, remains unverified. The published Loos et al. PRX appendix and equation locators are recorded in the reference index. Loos et al. and Schmiedl–Seifert are linked, not supplied as local PDFs. Large Zenodo trajectory archives are not included or claimed to have been independently reproduced.

## License

The manuscript and original code are distributed under **CC BY 4.0**. Third-party papers and data retain their original copyright, attribution, and license conditions. See [LICENSE.txt](LICENSE.txt).
