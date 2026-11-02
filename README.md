# Thermodynamic Accessibility and Hidden-State Structure in TCGS–SEQUENTION

Henry Arellano-Peña — manuscript version 4, 2 November 2026

## Research scope

This repository accompanies **Thermodynamic Accessibility and Hidden-State Structure in TCGS–SEQUENTION: Environmental Memory, Epistemic Information, and an Exact Finite-Time Work Law**.

The paper derives the exact finite-time accessible work for the specified passive, two-coordinate harmonic model with a free final trap centre. It distinguishes the full conditional free-energy bound, the hidden conditional resource, and the work accessible under declared observation and control restrictions. The implementation also supports the retrospective, no-branch-fit comparison with the source-data workbooks of Muruga and colleagues.

The complete distribution contains the unchanged 40-page version 4 manuscript (`main.pdf`), 27 numbered references, four figures, and seven tables. LaTeX sources and Overleaf archives are excluded. The mathematical results are conditional on the stated model; their TCGS–SEQUENTION interpretation is distinguished from the conventional operational explanation.

## Reproduction

Python 3.13.5 was used with the versions in `requirements.txt`.

```sh
python -m pip install -r requirements.txt
python tools/fetch_publisher_inputs.py
python code/reproduce.py
python code/label_resolution.py
```

The retrieval command accepts publisher files only if they match the recorded SHA-256 digests. Once the four original workbooks are present, reproduction does not require network access or random trajectory sampling.

- `code/accessibility.py`: exact Riccati solution, optimal mean path, constant-speed work, binary Gaussian conditioning, and independent discrete-control optimization.
- `code/reproduce.py`: original calculations, three figures, workbook summaries, and 70 numerical checks.
- `code/label_resolution.py`: complete ramp comparisons, all 24 label assignments, denominator reconstruction, conditional-mean magnitude inversion, optimal initial jumps, the fourth figure, and 35 additional checks.

Both suites pass: **70/70 and 35/35, totaling 105**. These are mathematical and numerical consistency checks, not 105 independent experimental tests. Units are nanometres, seconds, and zeptojoules. Physical parameters, threshold, delay, and protocols are taken from the source values identified in the manuscript; no branch-specific fitting is performed.

## Interpretation of the workbook comparison

The selected assignment is **B → 11, D → 00, F → 01, H → 10**. Original workbook headers and values are retained. This is a retrospective, no-branch-fit parameter-transfer comparison—not parameter-free physics, prospective validation, source-author confirmation, or independent reconstruction from raw experimental trajectories.

Conditional on transfer of stiffness and temperature, the column-D baseline residual corresponds to a +4.80320 nm conditional-mean magnitude offset; the corresponding F offset is +5.40195 nm. The inversion does not recover a sign. The denominator score separates D and H and identifies B/F as a pair; the full work curves resolve B versus F. The scores share model and data inputs and are not independent statistical confirmations. Instrument bandwidth is a possible, not established, explanation for the remaining residuals.

## Files and references

`references.json` and `REFERENCE_INDEX.html` identify all 27 numbered references. `REFERENCE_MANIFEST.sha256` identifies the exact reference PDFs and original workbooks. Other numerical outputs and verification reports are under `results/`.

The companion distribution supplies the compiled manuscript, four figure PDFs, twelve separate supporting-reference PDFs, and four original workbooks. A local path in the reference metadata identifies that exact distribution file; it does not imply that every binary has already been mirrored into this repository. The large Zenodo trajectory archives are not included or claimed to have been independently reproduced. Framework papers are distributed individually, not as compilation containers.

The exact 71-page *No Final Map* copy defines SR in §2.4, p.5. The public archive record is version 1, deposited 15 August 2026; equivalence to the attached undated reference copy remains unverified. Its exact identity is recorded in the checksum manifest. The Loos et al. PRX appendix and equation locators are provided in the reference index. Loos et al. and Schmiedl–Seifert are linked, not distributed as local PDFs.

## License

The manuscript and original code are distributed under **CC BY 4.0**. Third-party papers and data retain their own copyright, attribution, and license conditions. See `LICENSE.txt`.
