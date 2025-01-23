## Genome Assembly

### what is genome assembly
Genome assembly is the process of reconstructing the complete sequence of an organism's DNA from smaller fragments obtained through sequencing technologies. This process involves piecing together millions to billions of short or long DNA reads to form a continuous representation of the genome.

### What verkko fillet could do
There are many genome assembly methods available, but here we focus on using `verkko-fillet` to refine the results from the verkko assembler. `Verkko` is a graph-based tool that generates consensus sequences by accurately connecting sequence nodes. However, during this process, some nodes may be misconnected, or certain nodes may remain unconnected when they should be. `verkko-fillet` helps users identify such problematic nodes and edges, gather supporting information, and correct the connections to enable a subsequent verkko consensus run.

### Overall workflow
![overallworkflow](../../src/verkkofillet/data/test_giraffe/fig/giraffe_complete_verkko-genomeassembly.drawio.svg)

* Step 1: Run the assembler using the raw long-read dataset to generate the assembly.
* Step 2: Perform gap filling and graph fixing using `verkko-fillet`, then rerun the `verkko` consensus step to generate a new and improved consensus.
* Step 3: Identify and correct suspicious regions in the graph using error k-mers, trio markers, and coverage data. In this step, the focus is on fixing relatively larger regions or genomic blocks.
* Step 4: Use variant calling methods, such as DeepVariant or DeepPolisher, to detect and correct errors at a higher resolution.
