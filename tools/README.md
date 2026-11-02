# Publisher-input retrieval

`fetch_publisher_inputs.py` downloads the four original source-data workbooks from the publisher only when their SHA-256 digests match `REFERENCE_MANIFEST.sha256`.

```sh
python tools/fetch_publisher_inputs.py
```

The optional `--with-papers` flag also attempts to retrieve the published article, supplementary information, and public peer-review PDF. No alternate version is accepted. Existing files with a mismatched digest cause an error and are not overwritten.

The nine separately identified TCGS–SEQUENTION reference editions and the compiled manuscript are part of the companion distribution; the script does not replace those exact editions using newer or differently assembled public files.
