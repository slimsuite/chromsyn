# Frequently Asked Questions (and other tips)

## Interpreting ChromSyn Plots

### Q.	Red blocks indicate the reverse strand, however, if majority of the syntenic regions are on the opposite strand the chromosome is reversed and given the suffix "R". What is the difference between the two and are they contradicting or complementing each other?

One assembly is selected as the focal reference. Use `focus=X` to set this. Otherwise, the genome with the most "Forward" synteny with other genomes will be selected. (See the `#FOCUS` log message.) Scaffolds in other assemblies are then oriented to **maximise the forward synteny they share with their most syntenic scaffold hit** in another genome. If `orient=focus`, this will _always_ be their best hit in the focus genome. If `orient=auto` (the default), this orientation propagates out from the focus genome, but orientation (and scaffold order) is always set relative to the adjacent assembly (closest to the focus). If `orient=none`, no reversal takes place. Once the ordering and orientation of the scaffolds is set, the syntentic block are coloured blue or red _relative to the plotted orientation_. (You can see this from the crossover of the red blocks.) This means that a red syntentic block against one forward and one reversed (`*R`) scaffold will actually be a `+` strand synteny between the originals - see the synteny blocks in the Excel output for details if unsure.

### Q.	If telomeres (black circles) are significantly absent on scaffolds does this show errors in assembly? 

Telomeres are hard to assemble and will often be missing, especially for acrocentric chromosomes. Their absence merely indicates a lack of telomere-to-telomere assembly. Having telomeres at the ends of scaffolds demonstrates a good chromosome-level assembly, but their absence does not necessarily indicate errors. [Telociraptor](https://github.com/slimsuite/telociraptor) can help identify whether there has been mis-scaffolding at the ends of chromosomes that "hides" telomeres.

### Q. Is it possible to identify a centromere from the plots?

Currently, the only way to identify centromeres is to provide them through the features file. I intend to add specific centromere plotting in a future release, but identification will still need to be done with a different tool. AusARG have released one such approach, using TRASH, [here](https://github.com/kango2/ausarg?tab=readme-ov-file#centromeressh).
