# PTB
A density matrix (P) tight-binding (TB) method based on a polarized valence double-zeta basis set.
It is available for all elements and structures until Z = 86.

In this development version, before execution, two files have to be given:
The `.atompara` file containing all empirical parameters. If no specific location is given via a `-par` flag, it is assumed to be located in the `$HOME` directory (e.g., `~/.atompara`).
An individual location can be defined via the `-par <path of .atompara>` command.
The `.basis_vDZP` file contains the vDZP basis set in the correct format. Its file location can be defined via the `-bas` flag, otherwise it is assumed to be in 
`$HOME` (e.g., `~/.basis_vDZP`).

Molecular (total) charge can be incorporated via the presence of a `.CHRG` file. A similar procedure follows with the number of unpaired electrons (`.UHF`),
eventhough PTB is mainly developed for closed-shell systems.


## Building
You can use a statically-linked release binary (recommended), but you can also build it with the source code.
```
git clone 
cd source
```
You can the build the project via `make`.
