<table align="center">
  <tr>
    <td align="center">
      <img src="assets/istmed.png" width="400" alt="Sol Görsel" />
      <br />
      <sub><b>Sol Görsel Başlığı</b></sub>
    </td>
    <td width="40"></td> <!-- İki görsel arasındaki boşluk -->
    <td align="center">
      <img src="assets/durdagilab.png" width="400" alt="Sağ Görsel" />
      <br />
      <sub><b>Sağ Görsel Başlığı</b></sub>
    </td>
  </tr>
</table>

# DELDesign — reaction-constrained virtual DEL design

WP1 of the in-silico DNA-encoded library project. It enumerates only what
on-DNA chemistry and commercially available reagents could actually have made,
keeps the DNA attachment vector explicit throughout, and reports what a design
costs in validity and property compliance rather than only its nominal size.

```
pip install rdkit numpy
python tests/test_deldesign.py        # no test framework needed
python examples/design_demo.py        # compares five designs end to end
```

## Why this is not just combinatorial enumeration

The constraint that does the work is not the reaction SMARTS. It is the handle
bookkeeping: at every cycle the growing, DNA-bound intermediate must present
exactly one reactive handle, and the building block must both react now and
leave behind the handle the next cycle needs. A library definition that cannot
satisfy that chain is refused before anything is enumerated:

```python
>>> library.make_library("bad", "NCCCCCC*", [Cycle("c1", amide, amino_acids,
...                                                next_handle="amine")])
DesignError: cycle 'c1' carries an onward 'amine' while 'amide_coupling'
consumes the same handle: the reagent would self-couple. Declare the
protecting group (e.g. protecting_group='Fmoc') removed before the next cycle.
```

Three further decisions are made explicit rather than left implicit.

**Ambiguous reagents are rejected, not resolved.** A combination that could
give more than one product is dropped with the status `ambiguous_product`,
because a reagent whose on-DNA outcome a chemist could not predict does not
correspond to a single real well.

**The DNA tag is carried as a dummy atom and stripped only at scoring time.**
Profiling the on-bead conjugate makes every member look heavier and floppier
than the compound anyone will test — the amino-C6 headpiece alone contributes
six rotatable bonds. Declaring `linker_smarts` on the library removes the
headpiece for property calculations while the attachment vector stays in the
structure, which is what WP3 will need.

**Failures are part of the report.** A design that quietly loses a third of its
combinations at the last cycle and reports only the survivors looks far better
than it is, so `Profile` carries `validity`, the failure reasons, the cycle each
failure occurred in, and an `effective_size` — the nominal size scaled by the
fraction that is both valid and inside the property envelope.

## The pieces

| Module | What it holds |
|---|---|
| `reactions.py` | Seven on-DNA transformations, the handle SMARTS, per-reaction forbidden groups, and the reagent-level DNA hazard list |
| `blocks.py` | Catalogue curation: standardise, deduplicate on InChIKey, handle counting, hazard and property filters, with a per-reason rejection report |
| `library.py` | `Cycle` / `Library`, design validation, mixed-radix addressing, lazy construction |
| `profile.py` | Sampling-based property, scaffold, diversity and 3D shape profile; `compare()` for design review |
| `spec.py` | Library definitions as reviewable JSON |
| `cli.py` | `reactions`, `init`, `profile`, `sample`, `enumerate` |
| `data.py` | A small demonstration catalogue (174 reagents) so the pipeline runs without a catalogue licence |
| `examples/export_browser_data.py` | Enumerates the viewer's libraries in full |
| `examples/build_viewer.py` | Inlines data and RDKit into the standalone HTML |

## Addressing

A library is never materialised. It is defined by its headpiece, reaction
sequence and reagent sets, and any member is built on demand from its address —
the tuple of reagent indices the DNA tag would encode:

```python
lib.size                 # 99_179
lib.address(12345)       # (5, 2, 8)
lib.index_of((5, 2, 8))  # 12345
lib[12345].smiles        # built on demand
lib.sample(10_000, seed=7)     # reproducible, without replacement
```

So a 10^9-member design costs nothing to define, sampling is reproducible
across machines, and a hit found later can be traced back to the exact three
reagents that made it.

## Command line

```bash
python -m deldesign reactions                      # the reaction space
python -m deldesign init my_library.json           # a template spec
python -m deldesign init my_library.json --self-contained   # chemistry + reagents inline
python -m deldesign import my_design.json -o my_library.json   # from the HTML designer
python -m deldesign check   my_library.json        # curate and validate, build nothing
python -m deldesign profile my_library.json -n 20000 --json report.json
python -m deldesign sample  my_library.json -n 5000 --out sample.smi --address
python -m deldesign enumerate my_library.json --limit 200000 --out all.smi
```

A spec is a reviewable file; `builtin:` names the demonstration sets, anything
else is a path to `.smi`, `.csv` or `.sdf`:

```json
{
  "name": "PARP-focused trisynthon",
  "headpiece": "NCCCCCC*",
  "linker_smarts": "[#0]CCCCCC",
  "envelope": "lead-like",
  "cycles": [
    {"name": "c1", "reaction": "amide_coupling",
     "catalogue": "reagents/fmoc_amino_acids.sdf", "protecting_group": "Fmoc"},
    {"name": "c2", "reaction": "amide_coupling",
     "catalogue": "reagents/fmoc_amino_acids.sdf", "protecting_group": "Fmoc"},
    {"name": "c3", "reaction": "sulfonylation",
     "catalogue": "reagents/sulfonyl_chlorides.smi"}
  ]
}
```

The onward handle of each cycle is inferred from what the next cycle consumes,
so the spec cannot silently disagree with itself.

### Your own chemistry and your own reagents

Nothing in a design has to come from this package. A spec may define the
transformations it uses and carry its reagents inline, which makes it a single
file a collaborator can run and a reviewer can check:

```json
{
  "name": "my library",
  "headpiece": "NCCCCCC*",
  "reactions": {
    "my_acylation": {
      "smarts": "[N;H2,H1:1].[CX3:2](=[OX1:3])[OX2H1]>>[N:1][C:2]=[O:3]",
      "consumes": "amine", "reagent_handle": "carboxylic_acid",
      "regenerates": null, "conditions": "DMT-MM, pH 9.4 borate, rt"
    }
  },
  "cycles": [
    {"name": "c1", "reaction": "my_acylation",
     "reagents": ["OC(=O)c1ccccc1 benzoic", {"id": "nic", "smiles": "OC(=O)c1cccnc1"}]}
  ]
}
```

A reaction defined this way stays local to its spec: it does not enter the
package registry, so two designs may use the same name for different chemistry
without interfering. `python -m deldesign reactions` lists the handle names a
definition may use. If the group your reaction requires is also on the on-DNA
hazard list — a sulfonyl chloride is — the allowance is automatic.

`examples/design_your_own.py` is the same thing written as a plain Python
script, with the reactions, reagent lists and property limits at the top as
editable constants, and `--spec` to write the design out as JSON afterwards.

### Agreement with the browser designer

`DEL_designer.html` runs the same rules in the browser. On the worked
three-cycle design (22 × 18 × 14 = 5,544 members) the two front ends were
compared exhaustively, not by sampling: both enumerate 5,544 products, the
canonical on-DNA structure sets are identical member for member, and both
place 5,030 inside the lead-like envelope. A design exported from the page
therefore reproduces exactly here:

```bash
python -m deldesign import my_design.json -o my_library.json
python -m deldesign check my_library.json
```

## What the demo shows

Five designs built from the same cycle-1 reagents, profiled on equal terms:

```
design                            nominal   valid   envlp    effective  scaf/cmpd  meanTc  median MW
2-cycle / amide cap                 2,419 100.0%  89.5%        2,165       0.49   0.213     261.31
3-cycle / amide cap                99,179 100.0%  71.7%       71,078       0.82   0.212     381.44
3-cycle / sulfonamide cap          40,344 100.0%  50.7%       20,440       0.74   0.222     402.48
3-cycle / reductive-amination cap  36,982 100.0%  79.8%       29,523       0.82   0.221     359.47
3-cycle / SNAr cap                 18,491 100.0%  75.0%       13,868       0.71   0.225     350.93
```

The capping chemistry moves property compliance by nearly forty points at
constant reagent count, which is the kind of trade-off WP2 will optimise over
rather than choose by inspection.

## The DNA tag

A DEL is not a library of compounds but of compound-tag pairs, and the tag is
the only thing ever read. After a selection the compounds are gone; what
returns is sequencing reads, and the structures are recovered by decoding them.
`deldesign/encoding.py` designs the codons, builds the tags and reads them
back:

```bash
python -m deldesign encoding my_library.json --table --simulate 20000
python -m deldesign sample   my_library.json -n 1000 --dna --out tagged.smi
python -m deldesign decode   my_library.json --file reads.txt
```

```
[forward primer][codon 1][spacer][codon 2][spacer][codon 3][reverse primer]
```

Each cycle gets its own codon set. Four constraints are enforced rather than
assumed, and `encoding` reports all four: minimum Hamming distance between
codons, no homopolymer longer than two bases, GC held in a narrow band so
codons ligate at similar efficiency, and per-position base balance, without
which an Illumina run loses the contrast it needs to call clusters.

Substitutions are corrected, indels are only detected — by checking the
constant primer and spacer regions, which is what they are for. A read that
has slipped frame is discarded rather than assigned to whichever compound it
lands nearest.

### What the distance buys

Simulated on the worked three-cycle design, 40,000 reads per row, 15% of errors
as insertions or deletions:

```
min distance   error/base    correct   MISASSIGNED   rejected
           1        0.1%     97.485%        0.015%     2.500%
           1        0.5%     88.180%        0.022%    11.797%
           1        1.0%     77.683%        0.055%    22.262%
           3        0.1%     99.305%        0.000%     0.695%
           3        0.5%     96.315%        0.000%     3.685%
           3        1.0%     92.938%        0.007%     7.055%
```

The dominant cost of no error correction is read *loss*, not misassignment: at
0.5% per base, d = 1 recovers 88.2% of reads and d = 3 recovers 96.3%.
Misassignment is rarer by an order of magnitude but matters more per event, as
a lost read is merely lost while a misassigned one credits a count to a
compound that was never selected, and that is what sets the noise floor under a
weak enrichment.

### Codon capacity

Greedy selection under the filters above. A cycle needs one codon per building
block, and by default codons from different cycles are also held apart, which
draws the pool down once rather than once per cycle:

```
length     d=2     d=3     d=4     d=5
     5     120      25      10       4
     6     471      81      20       7
     7    1300     222      56      16
     8    5575     739     158      42
```

Distance 5 at length 7 gives only 16 codons, so a 22-member cycle cannot use
it; the error says so and names the fix. `distinct_across_cycles=False` relaxes
the cross-cycle constraint when a design is short of codons — position in the
read already separates the cycles, so this costs only robustness to a frame
slip.

### The same tag in both front ends

A tag is only an identifier if it is stable. The browser page assigns codons
with the same rules *and the same random number generator* as the package — a
small, fully specified xorshift rather than either language's built-in — and
the two were checked to produce byte-identical codon sets, layouts and tags for
the worked design. Python's generator cannot be reproduced in JavaScript, which
is why neither front end uses it.

## Synthetic accessibility

Ertl & Schuffenhauer's SA score (1 easy, 10 hard) is part of the descriptor set
and of the property envelope, bounded at 5 in `lead-like` and tunable per
design:

```json
{"envelope": "lead-like", "envelope_overrides": {"sa": [1.0, 4.0]}}
```

`LEAD_LIKE_NO_SA` lifts the bound entirely, for separating its effect from
everything else. `ReagentLimits(max_sa=...)` applies the same score to the
reagents instead of the products.

To see the score for every compound you generated, ask for the CSV:

```bash
python -m deldesign sample    my_library.json -n 2000 --csv --dna --out sample.csv
python -m deldesign enumerate my_library.json --csv --out all.csv
```

Each row carries the off-DNA and on-DNA SMILES, the building block used in each
cycle, the DNA tag, every descriptor including `sa`, and whether the compound
falls inside the envelope. `--props` adds the same columns to the plain SMILES
output instead.

**It is nearly inert on DEL products, and that is a property of DELs, not a
bug.** Every member is built from commercial reagents by two or three simple
bond formations, which is exactly the kind of molecule the score was fitted to
call easy. Measured on 1,200-member samples of five designs sharing one
cycle-1 reagent set:

```
design                              SA median   SA max   > 5   > 4   > 3
2-cycle / amide cap                      2.22     4.41     0     7   4.1%
3-cycle / amide cap                      2.66     4.46     0    16  19.1%
3-cycle / sulfonamide cap                2.78     3.86     0     0  27.9%
3-cycle / reductive-amination cap        2.70     3.65     0     0  23.1%
3-cycle / SNAr cap                       2.78     3.67     0     0  29.6%
```

Not one compound in six thousand exceeds 5. On the full 5,544-member
enumeration the cut at 5 removes zero members and the lead-like fraction stays
at 90.7%.

The score is not useless here, but what it is useful *for* is worth being
precise about. It says nothing about whether the library can be made — the
reaction and handle constraints in this package are what decide that. What it
estimates is the cost of off-DNA resynthesis, the step without which a DEL hit
cannot be confirmed or optimised. That makes it a statement about hit
follow-up.

A cut tight enough to bind works against shape diversity. On the same
enumeration:

```
SA cut   removed        mean fsp3   mean aromatic rings   mean stereocentres
none        0  ( 0.0%)      0.386                  1.89                 0.41
 3.0      750  (13.5%)      0.348                  2.02                 0.38
 2.8    1,801  (32.5%)      0.323                  2.11                 0.33
 2.6    3,177  (57.3%)      0.294                  2.19                 0.24
```

SA correlates positively with sp3 fraction (r = +0.57) and negatively with
aromatic ring count (r = −0.44): it rewards flat molecules. Tightening it
therefore undoes what mixing the backbone chemistry was for. Use it on the
reagents, where a real catalogue does contain awkward entries, rather than on
the products.

### In the browser too

The page computes the same score, under every structure. RDKit's WebAssembly
build exposes only folded fingerprints while this score needs unfolded Morgan
identifiers with counts, so the fingerprint is reimplemented in JavaScript
from the molecule graph — which is only worth doing if it is exact, so it was
checked against RDKit rather than assumed: **identical scores to eight decimal
places on all 5,544 members** of the worked library, and on the published
reference values (aspirin 1.58, nicotine 2.50, estradiol 3.59, morphine 5.27).
Finding that took four separate bugs, each of which produced plausible but
wrong numbers: the per-layer hash combines each neighbour *pair* as one value
rather than its two members; potential stereocentres must be counted on a
stereochemistry-stripped copy, because a trans-1,4-disubstituted cyclohexane
carries two CIP labels and no chiral centres; bridgehead atoms arise only when
two rings share three or more atoms; and an environment comparator that
indexed the wrong array scrambled the duplicate elimination.

## Scaling to a thousand reagents per cycle

`deldesign/catalogue.py` generates large reagent sets by decorating cores with
one or two substituents, and every handle can supply well over a thousand:

```bash
python -m deldesign catalogue                    # what is available
python -m deldesign catalogue aldehyde -n 1000 --out aldehydes.smi
```

```json
{"name": "c1", "reaction": "amide_coupling", "catalogue": "generated:amino_acid:1000"}
```

The same sets are embedded in `DEL_designer.html`: each cycle has a
**Fill with building blocks** button offering 20, 200, 600 or 1,400 reagents,
and the page curates 3 × 1,400 of them in about nine seconds. Because a
thousand building blocks per cycle cannot be encoded by a 7-mer, the codon
length there defaults to *chosen automatically*: it finds the shortest codon
that can supply the design, and only if none up to ten nucleotides can will it
give up holding different cycles' codons apart, saying so when it does. For
1,177 × 1,382 × 1,337 it settles on 10-mers and a 56 nt tag.

**These are virtual reagents, not a procurement list.** They are chemically
reasonable structures produced by combinatorial substitution; some fraction of
any set is not purchasable at any price. They exist to size a design and to
test that a pipeline survives a realistic number of building blocks. A real
design replaces them with a supplier catalogue, through the same curation.

Two consequences of asking for a thousand per cycle are worth seeing before
committing to it.

**The codon has to grow.** A thousand building blocks needs a thousand codons
per cycle, and 8-mers yield only 739 under the GC, homopolymer and distance
filters. Nine is the shortest that works, which takes the tag from 47 to 53
nucleotides; holding codons apart *across* cycles as well needs ten.

**The returns are strongly sublinear.** Reagent sets reach a thousand by
getting larger and more decorated, and product size is close to the sum of
reagent sizes, so the products leave lead-like space. Measured on 900-member
samples of the same three-reaction design, with reagent limits already
tightened to MW 260, cLogP 3, 5 rotatable bonds and 19 heavy atoms:

```
per cycle         nominal   valid  lead-like      effective  med MW  med SA
       25          15,625   91.7%      33.2%          4,756  447.56    2.89
      100       1,000,000   93.9%      15.7%        147,777  474.56    3.04
      400      64,000,000   95.2%       8.6%      5,262,222  496.03    3.78
     1000   1,000,000,000   89.8%       3.2%     28,888,888  518.68    3.85
```

The browser reproduces this independently: filling all three cycles with 1,400
reagents gives 1,177 × 1,382 × 1,337 = 2,174,782,918 members, of which 90% are
chemically valid and **1% are lead-like** — an effective size of about 14
million. Going from 25 to 1,000 per cycle multiplies the nominal size by
64,000 and the effective size by only 6,075. The binding constraint is the property envelope,
not the reagent supply, and the way to spend a larger reagent budget well is on
*smaller* reagents rather than more decorated ones. Note also that the SA
score only begins to bite at this scale: its median moves from 2.89 to 3.85 and
the upper tail finally reaches the cut at 5.

## Measuring a design against compounds that are known to work

`examples/compare_to_known.py` compares a library with actives, and ships with
the IL17A series from J. Med. Chem. 2024, 67, 6456-6494 (Ramos et al.), found
by DEL screening at X-Chem/AbbVie. That campaign is the right yardstick here
because its libraries are close cousins of what this package builds: 746
Fmoc-amino acids acylated onto an amine-terminated linker then 4,920 reagents
installed by acylation, urea formation or sulfonamide formation; and two
cycles of Fmoc-amino acids capped with 2,462 carboxylic acids.

Against the worked three-reaction design, the answer is a clean negative:

```
hit   MW     ours (amide/RA/sulfonylation)   paper topology, our reagents
 1   354.2              0.338                           0.492
 5   320.2              0.240                           0.404
 6   332.2              0.228                           0.268
 7   493.4              0.196                           0.247
 8   628.6              0.193                           0.237
```

Maximum Tanimoto 0.19-0.34 (Morgan r2, 2048 bit) means the library contains
nothing of this kind. Two separate causes, and they are worth separating.

**The reaction scheme.** Both published hit series are capped peptides — pure
acylation chemistry. Replacing reductive amination and sulfonylation with two
more amide couplings, *with the same reagents*, moves the closest match from
0.338 to 0.492 for compound 1. Mixing the backbone chemistry widens shape
diversity, which is a real gain, and it costs access to this hit class.

**The reagents, which matter more.** Six of the seven building blocks the paper
names as enriched are absent from the generated catalogue, including the
4,5-dichloro-1H-pyrrole-2-carboxamide that every member of the series carries.
A pharmacophore that is not in the reagent list cannot be enumerated at any
library size. Given the paper's own reagents and its own all-acylation
topology, this engine reaches the hits from **27 compounds**:

```
hit 1  Tc 0.603    hit 5  Tc 0.435    hit 6  Tc 0.750
hit 7  Tc 0.573    hit 8  Tc 0.816
```

Curation accepts all of those reagents in their correct cycle roles, so the
filters are not what excluded them — the catalogue is.

**The envelope would also have discarded the better hits.** Compounds 1, 5 and
6 sit inside the lead-like window, but 7 fails on molecular weight (493) and 8
fails on weight, cLogP and heavy atoms (629 Da, cLogP 5.3). Compound 8 came
from the three-cycle library and is the more advanced starting point. A
property window tuned for lead-likeness removes the better protein-protein
interaction hits, which is the known tension between lead-like filtering and
PPI targets rather than a flaw in either.

## The DNA attachment point in docking

During a selection the compound is never free: it hangs off a duplex roughly
20 Å across. Docking the off-DNA compound forgets this, and the usual
consequence is a pose that could not have existed — the attachment atom buried
at the bottom of the pocket with the tag pointing into the protein. The
failure is silent: the pose scores well and nothing in the output says the
molecule was tethered.

`deldesign/vector.py` keeps that information and uses it, in four steps.

**Dock the species you will resynthesise, not the bare fragment.** The
compound a chemist makes off-DNA carries a small group where the linker was —
the IL17A paper above used a methylamide on compound 1 and an N-phenylacetamide
on compound 8. `--cap amide` or `--cap methylamide` builds that molecule rather
than a truncated fragment, so the sterics and the hydrogen bonding at the
attachment position are the ones that will be tested.

**Record which atom held the DNA.** It is not the dummy and not the first
linker atom — both vanish on truncation — but the first atom that survives,
usually the headpiece amide nitrogen. It is written into the SD file as
`del_attachment_atom` so it survives the round trip through a docking program.

```bash
python -m deldesign dock-prep my_library.json -n 2000 \
       --cap methylamide --dna --out for_docking.sdf
```

**Constrain the docking itself where the software allows it.** Published DEL
docking work places a positional constraint at the solvent-exposed entrance of
the site so poses must present the exit vector outward. Glide takes positional
and excluded-volume constraints directly; MOE takes pharmacophore exclusion
volumes; with AutoDock Vina or smina there is no native equivalent, so the
filter below does the work after the fact.

**Filter the poses on where the tag could go.**

```bash
python -m deldesign vector-check --poses docked.sdf --receptor receptor.pdb \
       --accept open --out keepers.sdf
```

For each pose this walks outwards along the exit vector in 1 Å steps and
reports how far it gets before a receptor heavy atom comes within 2.8 Å, plus
how many receptor atoms sit inside a 35° cone around it. `clearance` is
reported in ångströms rather than as an arbitrary score, so it can be
inspected and argued with: `open` means the vector runs the full 12 Å, `tight`
means at least 6, `blocked` means the tag has nowhere to go. Twelve ångströms
is about the span of the linker alone, so a vector that is not open at 12 Å
should be treated as unavailable rather than merely tight.

Two limits worth stating plainly. The thresholds are geometric defaults, not a
validated protocol — the published accessibility scores for DEL docking are
not standardised, so rank and inspect rather than cutting blindly. And the DNA
is not always innocent: it can contribute to binding, and attachment vectors
are documented to interfere with ternary complex formation in degrader work,
so a pose that passes this test is not proved, only not excluded.

## Developability: structural alerts and property rules

```bash
python -m deldesign admet my_library.json -n 5000
python -m deldesign sample my_library.json -n 2000 --csv --admet --out s.csv
```

The page has the same filters: **No reactive or liability groups**, and a
*must pass* menu for Lipinski, Veber, Egan, Pfizer 3/75, GSK 4/400 or all of
them. Each structure shows its flags and an estimated logS, and the CSV gains
`logS`, `n_alerts`, `alerts` and `rules_failed`. Page and package were
cross-checked on 1,857 compounds: **no disagreement on any alert, logS value
or rule verdict**.

### Trained models, as an optional extra

An earlier version of this file said organ toxicity could not be predicted
here. That was too broad, and it is corrected: trained models for most of
these endpoints exist, ship their own weights and run offline.

```bash
pip install admet_ai                      # pulls in PyTorch
python -m deldesign predict my_library.json -n 2000 --check --out preds.csv
```

This uses ADMET-AI (Swanson et al., *Bioinformatics* 2024), Chemprop-RDKit
ensembles trained on Therapeutics Data Commons. It covers hERG, DILI, Ames,
carcinogenicity, ClinTox, skin sensitisation, LD50, CYP3A4/2D6/2C9 inhibition,
microsomal clearance, half-life, plasma protein binding, volume of
distribution, Caco-2, intestinal absorption, oral bioavailability, P-gp, BBB,
solubility and logD — at about 75 compounds a second, so a 2,000-member sample
takes half a minute and a million would take four hours.

`--check` runs reference compounds first, because a model that silently
changes version is worse than no model:

```
terfenadine (hERG blocker)    hERG 0.97     atenolol (low risk)  hERG 0.14
astemizole  (hERG blocker)    hERG 0.99     aspirin  (low risk)  hERG 0.02
ketoconazole (hepatotoxic)    DILI 0.93     metformin            DILI 0.00
benzo[a]pyrene (mutagen)      AMES 0.90     caffeine             AMES 0.11
```

hERG and Ames separate cleanly, and hERG catches sertindole (0.77), which the
lipophilicity heuristic misses.

**Three caveats that the package states rather than hides.** The **DILI model
is weak**: it scores caffeine 0.93, aspirin 0.67 and benzo[a]pyrene 0.94, so a
high DILI score is close to uninformative on its own and a low one is the more
meaningful signal. **The regressions extrapolate**: on 1,200 members of the
worked library, predicted half-life was *negative* for 77% of them, plasma
protein binding came back above 100%, and volume of distribution landed in the
5th percentile of approved drugs. Rank within a library with these; do not read
them as absolute values. And **percentiles against approved drugs are the
number to read first** — a library with median DILI 0.60 sits at the 61st
percentile, which is where a great many ordinary molecules sit.

**Still not covered by any model here: nephrotoxicity and cytotoxicity.**
ClinTox and LD50 are the nearest stand-ins and are different endpoints.

### The page itself contains no organ-toxicity prediction

No hepatotoxicity model, no nephrotoxicity model, no cytotoxicity model, no
trained cardiotoxicity model. Not one of the 301 alert names in the
catalogues used refers to an organ or an endpoint — they are names of chemical
groups (`Michael_acceptor_1`, `N-nitroso`, `Polycyclic_aromatic_hydrocarbon`).
A Chemprop ensemble cannot run in a browser, so the page keeps the heuristics
and the structural alerts, and says so. `admet.NOT_PREDICTED` lists what the
page does not do, so the absence cannot be inferred from silence:

```
hepatotoxicity (DILI) · nephrotoxicity · cytotoxicity · cardiotoxicity/hERG as
a trained model · mutagenicity (Ames) · carcinogenicity · clearance, half-life,
bioavailability · CYP inhibition or induction · plasma protein binding
```

Each needs a model trained on measured outcomes — DILIrank records, Tox21 cell
viability, hERG patch-clamp currents — which is a different kind of object from
a substructure match.

Two toxicity-related signals *are* present, both correlations drawn from sets
of marketed drugs rather than predictions:

* **Pfizer 3/75** (Hughes et al., Bioorg Med Chem Lett 2008): cLogP above 3
  together with TPSA below 75 went with markedly more preclinical toxicity.
* **hERG structural risk**: a protonatable basic amine plus lipophilicity
  (cLogP > 3.7) is the recurring motif among potent blockers. Against ten
  reference drugs it flagged four of five well-known blockers — terfenadine,
  astemizole, haloperidol, verapamil — and none of five low-risk compounds —
  atenolol, metformin, aspirin, paracetamol, ranitidine. It missed sertindole,
  whose Crippen cLogP of 2.5 badly underestimates the measured value. A flag is
  a reason to measure; no flag is no information.

Separately, `admet.mechanisms()` groups the alerts by the bioactivation route
they are associated with — nitroreduction, aromatic-amine N-hydroxylation, ring
epoxidation, Michael adduction, quinone redox cycling, direct alkylation. That
adds no prediction; it says why a group is on the list, which is the difference
between a flag you can act on and one you cannot.

**This is not ADME prediction, and nothing here claims to be.** The IL17A
paper makes the case better than any argument: compounds 1 and 8 share a
scaffold and a pharmacophore, and differ by roughly a thousandfold in unbound
human microsomal clearance (under 2.4 against 2100 L/h/kg) and three
hundredfold in PAMPA permeability (11 against 0.032 x 10^-6 cm/s). No 2D
filter separates those correctly. What is decidable before synthesis is
narrower: whether a compound carries a group with a known liability, and
whether it sits where marketed drugs mostly sit.

### Screen the products, not the reagents

This is the opposite of the advice for synthetic accessibility, and the reason
is specific to DELs: **the commonest alerts on DEL reagents are the reactive
handles themselves.** Across a 643-reagent set for one three-reaction design,
the alerts were 440 aldehydes, 220 sulfonyl halides and 220 halogen-heteroatom
bonds — every one consumed by the chemistry and absent from every product.
Filtering reagents naively would delete most of the catalogue for carrying the
groups the reactions need. The alerts that actually survive into products are
different: nitro groups, dialkyl anilines, alkyl halides, 2-halopyridines.

Measured on 1,500 products of that design: 80.1% carry no alert, 16.9% carry
only alerts inherited from a building block, and 3.0% carry one the coupling
itself created. So reagent screening still helps, but only after discounting
the handle the reaction requires — `admet.reagent_alerts` does exactly that,
and leaves a genuine alert elsewhere on the same reagent in place.

### What runs where

The package runs RDKit's full catalogues — BRENK, NIH, PAINS A/B/C, ZINC and
the ChEMBL sets, hundreds of patterns — selectable with `--catalogues`. The
page runs a named list of fourteen reactive and liability groups, because
RDKit does not expose its catalogue patterns for export and transcribing
hundreds of SMARTS from a paper would be a quiet source of error. Each of the
fourteen was checked to fire exactly when RDKit's entry of that name fires,
across hundreds of molecules, and they are listed in `admet.VERIFIED_ALERTS`
rather than hidden. The page says plainly that it checks fourteen and the
package checks far more.

Rules are reported as the specific criterion broken, not as a pass/fail: a
compound does not "fail Lipinski", it exceeds a weight or a cLogP, and which
one decides whether it matters.

## Limits worth stating

A reaction may declare `regenerates`, the handle it leaves at the reaction site:
a reductive amination turns the amine it consumes into a secondary amine, so the
chain continues without the reagent carrying an amine of its own and without a
protecting group. Encoding that is what allows consecutive cycles to use
different chemistries. The reaction set is representative, not exhaustive;
adding a transformation is one dataclass entry. Yields are not modelled, so `validity` is a structural
measure of whether a product is well defined, not a prediction that the well
will contain it. The demonstration catalogue is far too small to design a real
library and exists only so the pipeline runs. Protecting-group handling is
recorded rather than simulated: the enumerated structure is the deprotected
product, which is what the DNA tag encodes.

## The browser front ends

`DEL_designer.html` is the interactive one: you pick or write the reaction
(SMARTS, the handle it consumes, the group the reagent must carry, optionally
the handle it regenerates), paste your own building blocks, and the page
curates with a per-reason rejection report, validates the handle chain,
enumerates, profiles, draws the structures and downloads CSV, SMILES or the
design itself. A downloaded design loads back into the page, or into this
package with `python -m deldesign import`. Rebuild it with
`python viewer/build_designer.py`.

`DEL_viewer.html` is the fixed one. It carries four fully enumerated libraries — every structure, every
descriptor — plus a compiled copy of RDKit, so it opens by double-clicking and
works with no internet connection and nothing installed. Pick a library, read
how it was made, filter with sliders, look at the structures, and download the
selection as CSV or SMILES. `START_HERE.txt` is the one-page guide for a
non-specialist.

Rebuild it after changing the reagents or designs:

```bash
python examples/export_browser_data.py    # enumerate, write browser_data.json
python examples/build_viewer.py           # inline everything into the HTML
```

The RDKit WebAssembly build is read from `/tmp/rdkit_inline.js` and
`/tmp/rdkit_wasm.b64`; regenerate them with
`npm install @rdkit/rdkit` and base64-encode `dist/RDKit_minimal.wasm`.
Note that this build ignores `Module.wasmBinary`, so the viewer supplies the
module through the `instantiateWasm` hook instead — without that it tries to
fetch the `.wasm` over the network, which a browser refuses to do from a local
file.

## Next

WP4 — hierarchical virtual screening of an enumerated library against a target,
with an active-learning funnel and a pluggable docking backend — consumes the
output of this package directly.

---

## GitHub release and development

This project is an **experimental research prototype**, not an experimentally
validated DEL synthesis protocol or an affinity predictor. Reaction SMARTS,
DNA compatibility, and individual products require independent chemical and
laboratory validation.

### Install

```bash
python -m venv .venv
source .venv/bin/activate       # Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install -e .
python -m deldesign --help
```

The project requires Python 3.10+ and the dependencies in `pyproject.toml`.
Optional trained ADMET models require `pip install -e '.[models]'`.

### Tests

```bash
python -m pip install pytest
python -m pytest -q tests
python -m compileall -q deldesign examples viewer
```

GitHub Actions runs Python tests on Python 3.10, 3.11 and 3.12.
The `viewer/test_designer.js` and `viewer/test_edge.js` files are **optional**
Playwright browser checks. They require a separately installed Node.js
Playwright environment and browser binaries, and are not part of the Python CI.

### Standalone offline browser interface

Open `DEL_designer.html` in a browser. The generated HTML embeds RDKit WebAssembly
and needs no application server. To rebuild it reproducibly from the checked-in
HTML and the editable HTML shell:

```bash
python viewer/build_designer.py --output /tmp/DEL_designer.html
```

To use your own locally obtained RDKit.js build instead of extracting the
bundled runtime, provide `--assets-dir /path/to/rdkit-js` containing
`RDKit_minimal.js` and `RDKit_minimal.wasm`. The default mode relies on the
existing checked-in `DEL_designer.html` being present; preserve it in the repo.

### Licensing

The repository owner maintains the **MIT LICENSE** in the GitHub repository.
Before publishing, make sure that file is included at the repository root.
RDKit and any bundled third-party components retain their own licenses. Review
and preserve the corresponding third-party copyright and license notices.

### Releasing

Treat initial releases as pre-release (`v0.1.0`). Validate the browser interface,
reaction outcomes and representative examples before claiming production readiness.

### GitHub Pages

The repository includes `index.html`, which redirects to the standalone
`DEL_designer.html`. In repository **Settings > Pages**, set **Deploy from a
branch**, choose `main` and `/ (root)`. GitHub Pages will serve the app from
`https://OWNER.github.io/REPO/` after publishing is enabled. Test it in a real
browser after deployment; the included HTML uses a large embedded WASM module.

### Test policy

The default CI executes the 14 short Python tests; the large reagent-generation
`test_catalogue` is excluded from push/PR CI because it takes substantially
longer. To run the entire 15-test suite, go to GitHub **Actions > Full chemical
catalogue tests > Run workflow**, or run `python -m pytest -q tests` locally.
A complete successful run of the large catalogue test has **not** been confirmed
in the present packaging review.
