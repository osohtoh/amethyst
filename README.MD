# amethyst - an analogue generator
---

It takes core molecule and list of R-groups and returns all possible combinations. For now documentation is available in the code via docstrings with more in-depth one being planned for the future.

Requires rdkit for core functionality and click for CLI access.

## CLI usage

`generate` command

```
Usage: amethyst generate [OPTIONS] CORE_SMI R_FILE

  Generate combinations of provided CORE_SMI and R_FILE.

Options:
  -d, --delimiter TEXT        Symbol separating R-groups.  [default: ,]
  -e, --enantiomers           Option to generate possible enantiomers.
  -o, --output-file FILENAME  File output path. Provided filepath must be
                              empty.
  --debug                     Enable debug logs. Funnily enough currently
                              broken so it always debug logs.
  --help                      Show this message and exit.```
