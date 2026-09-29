# English release verification

The release retains the general geological modelling workflow and documents Munich as a reference example. Chinese-language README files, internal editorial comments and the D2.2 review draft are excluded. The English D2.1 reference remains available.

Only comment text was translated in the numerical utility and legacy workflow. Before writing the translated files, every Python token other than comments was compared with the original source for executable modules; exact non-comment lines were compared for the syntactically incomplete legacy reference. Original hashes, release hashes and the unchanged executable-token digests are recorded in `core_source_localization.json`. Historical simulation manifests keep their original numerical-source hashes.

The regression suite includes a check against the original executable-token digest. Machine-local workspace path prefixes were replaced with `<LOCAL_WORKSPACE>` in selected JSON provenance records and the QA report source-file column; these are historical metadata, not runtime paths. Numerical arrays, observations, configuration values, figures and VTK data remain unchanged.

Final regression, reference-evidence, language, local-link and file-integrity results are recorded in `english_release_checks.json`. Complete production and cross-validation simulations were not rerun as part of language cleanup. The supplied D2.1 document is retained without content or layout changes.
