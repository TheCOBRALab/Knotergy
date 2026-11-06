# RNAstructure parameter sets

The nearest-neighbor parameters efn2 uses, from RNAstructure 6.6 (April 2, 2026), with one folder
per set. The files are unchanged copies of `RNAstructure/data_tables` and keep RNAstructure's file
names (`<prefix>.<table>.dg` for free energies at 37 °C, `.dh` for enthalpies), so for
pseudoknot-free structures a folder can also be efn2's `DATAPATH`, e.g.
`DATAPATH=params/RNAstructure/m1Ψ efn2 --alphabet 1 ...`.

| Folder        | efn2 `--alphabet` | Parameters                                     | Extra symbols          | Other temperatures |
| ------------- | ----------------- | ---------------------------------------------- | ---------------------- | ------------------ |
| `turner2004/` | `rna`             | RNA, Turner 2004 rules (default)               |                        | yes                |
| `dna/`        | `dna`             | DNA                                            |                        | no (1)             |
| `dmem/`       | `dmem`            | RNA in *in vivo*-like conditions               |                        | no                 |
| `m6A/`        | `m6A`             | RNA + N6-methyladenosine                       | `M`, `6`, `m` = m6A    | no                 |
| `Ψ/`          | `PSU`             | RNA + pseudouridine                            | `P`, `p` = Ψ           | no                 |
| `m1Ψ/`        | `1`               | RNA + N1-methyl-pseudouridine                  | `1`, `M`, `m` = m1Ψ    | no                 |
| `Ψ_full/`     | `FullP`           | RNA in which every U is Ψ                      | `U`, `T`, `P` = Ψ      | no                 |
| `m1Ψ_full/`   | `Full1`           | RNA in which every U is m1Ψ                    | `U`, `T`, `1`, `M` = m1Ψ | no               |
| `pseudoknot/` |                   | `pseudconst.dat`: efn2's pseudoknot penalties  |                        |                    |

Every set has the bases `A`, `C`, `G`, `U` (= `T`); lower-case letters are forced unpaired, `X`
and `N` neither pair nor stack, and `I` is the intermolecular linker. Each folder is a complete
set (like `efn2 --alphabet m6A`); Knotergy can also add the m6A, Ψ and m1Ψ sets to a base set
such as `turner2004`, several at once. The `_full` sets only make sense as the base set.

(1) Only free energies at 37 °C are available for the sets without `.dh` files. `dna` has them
for every table except `helix_ends`, so RNAstructure 6.6 cannot rescale it either.

## Sources

According to RNAstructure's documentation (`manual/Text/alphabets.html`,
`manual/Overview/Parameters.html`, `manual/Overview/Updates.html`):

- `turner2004`: the Turner 2004 rules as defined by Mathews et al. (2004) PNAS 101:7287-7292,
  with multibranch loop parameters revised to experimentally determined values (Lu et al. 2009).
  Enthalpies: Lu, Turner & Mathews (2006) Nucleic Acids Res. 34:4912-4924. In RNAstructure 6.6
  the terminal-pair penalties also depend on the penultimate pair (`rna.helix_ends`).
- `dna`: fit to available experimental data.
- `dmem`: fit to optical melting experiments in advanced DMEM buffer, based on the Turner 2004
  rules (added in RNAstructure 6.6).
- `m6A`: Kierzek et al. (2022) Nat. Commun. 13:1271; Szabat et al. (2022) J. Mol. Biol.
  434:167632.
- `Ψ`, `m1Ψ`, `Ψ_full`, `m1Ψ_full`: added in RNAstructure 6.6; the documentation cites no
  publication for them.
- `pseudconst.dat`: the constants efn2 reads when a structure has a pseudoknot (the file
  describes them as ShapeKnots' pseudoknot penalty).

RNAstructure is distributed under the GNU General Public License, version 2.

Not copied: `rna.cov.dg/.dh` and `rna.param_map.dg` (10 MB; efn2 only uses them for the ± error
estimate it prints for `rna`), `dynalignmiscloop` tables (Dynalign), `msa` (alignments), and the
`*-cuda` legacy tables.

## Adding a set

Make a folder holding `<prefix>.specification.dat` and the 20 tables `<prefix>.<table>.dg`
(`.dh` files are only needed at other temperatures). The prefix is optional (`stack.dg`,
`specification.dat`), and the folder's name is the set's name in Knotergy:

```c++
// Turner 2004 + m6A + m1Ψ at 37 °C; RNAstructure's names work too ("rna", {"m6A", "1"}).
const rnastructure_param p = RNAstructureParams::load("turner2004", {"m6A", "m1Ψ"});
```
