# Theseus-lib

**TABLE OF CONTENTS:**

* [Introduction](#introduction)
    * [What is Theseus?](#what_is_theseus)
    * [Getting started](#getting_started)
* [Using Theseus](#using_theseus)
    * [Multiple Sequence Alignment (MSA)](#msa_usage)
    * [Aligning a sequence to a graph](#alignment_usage)
    * [Creating a graph](#graph_creation)
    * [Penalties and dual-affine gaps](#penalties_usage)
* [Tools](#theseus_tools)
* [Heuristics](#theseus_heuristics)
* [Datasets](#theseus_datasets)
* [Reporting Bugs and Feature Request](#theseus_bugs)
* [License](#theseus_licence)
* [Authors](#theseus_authors)
* [Citation](#theseus_cite)

## <a name="introduction"></a> 1. Introduction

### <a name="what_is_theseus"></a> 1.1. What is Theseus?
Theseus is a fast, optimal and affine-gap Sequence-to-Graph aligner. It leverages the expected high similarity in the aligned data to accelerate computation and reduce the search space compared to other alternatives. Theseus is a general purpose aligner, providing two main functionalities:
1. **Multiple Sequence Alignment (MSA):** Theseus performs MSA of a set of N sequences using the Partial Order Alignment (POA) approach. That is, it progressively builds a partial order graph representing an MSA of a set of given sequences, adding one more sequence to the graph at each iteration.
2. **Alignment to a reference graph:** Alternatively, Theseus can align one or several sequences to a reference graph. The user has to provide an initial position for the alignment to start.

<p align = "center">
<img src = "img/Theseus_green.png" width="300px">
</p>

Theseus extends the proposal from the Wavefront Alignment algorithm (WFA), originally devised for pairwise sequence alignment, to the context of sequence-to-graph alignment.

### <a name="getting_started"></a> 1.2. Getting started
Git clone and compile the library, tools, and examples (by default, use `cmake` for the library and benchmark build).

```
git clone https://github.com/albertjimenezbl/theseus-lib
cd theseus-lib
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```


## <a name="using_theseus"></a> 2. Using Theseus

### <a name="consensus_usage"></a> 2.1. Multiple Sequence Alignment

This example illustrates how to use Theseus as a Multiple Sequence Aligner. First, include the msa and general headers:
```
#include "theseus/alignment.h"
#include "theseus/heuristics.h"
#include "theseus/penalties.h"
#include "theseus/theseus_msa_aligner.h"
```

Then, create and configure a MSA aligner object. This object is defined by three parameters: a set of penalties, a heuristics object, and an initial sequence behaving as the starting point of the alignment. An example on how to set such parameters and create an MSA object is found in the next code snippet:
```
theseus::Penalties penalties(match, mismatch, gap_open, gap_extend);
theseus::Heuristics heuristics();
theseus::TheseusMSA aligner(penalties, heuristics, initial_sequence, seq_weight);
```

Once this is done, we can start adding sequences to our POA graph using the align functionality, that returns an Alignment object with CIGAR, path, and score information:
```
theseus::Alignment alignment_object = aligner.align(sequence, weight, lag_pruning_active);
```

Each time a new sequence is added to the POA graph, the graph is updated with the newly found variation (all the insertions, deletions and mismatches of the resulting alignment object).

Finally, we can output the result in five different formats: a graph in .gfa format, a multiple sequence alignment, a consensus sequence, a consensus sequence based on the majority voting algorithm, and a graph in .dot format:
```
// Output as a Multiple Sequence Alignment
aligner.print_as_msa(output_file);

// Output the compacted POA graph in .gfa format
aligner.print_as_gfa(output_file);

// Find the consensus sequence of the alignment (Heaviest bundle algorithm)
sequence = aligner.heaviest_bundle_consensus();

// Find the consensus sequence of the alignment (Majority voting algorithm)
std::vector<int> consensus_weights;
std::string consensus_sequence, consensus_sequence_gapped;
aligner.majority_voting_consensus(consensus_weights, consensus_sequence, consensus_sequence_gapped);

// Output the compacted POA graph in .dot format
aligner.print_as_dot(output_file);
```


### <a name="alignment_usage"></a> 2.2. Aligning a sequence to a graph

This example illustrates how to use Theseus as a general Sequence-to-Graph Aligner. First, include the alignment and general headers:
```
#include "theseus/alignment.h"
#include "theseus/graph.h"
#include "theseus/penalties.h"
#include "theseus/heuristics.h"
#include "theseus/theseus_aligner.h"
```

Then, create and configure an aligner object. This object is defined by three parameters: a set of penalties, a heuristics object, and a reference graph in Theseus' internal graph format. An example on how to set such parameters and create an aligner is found in the next code snippet:
```
theseus::Penalties penalties(match, mismatch, gap_open, gap_extend);
theseus::Heuristics heuristics();
theseus::TheseusAligner aligner(penalties, heuristics, std::move(graph));
```

Once this is done, we can start aligning sequences to the reference graph using the align functionality. **Importantly**, a call to the align function consists of five arguments: the **sequence** to be aligned, the **starting vertex** for the alignment, the **starting offset** in that starting vertex, and two boolean variables indicating the **usage of the drop heuristic**, and the **usage of the lag pruning heuristic**. The result of the alignment is an Alignment object with CIGAR, path and score information:
```
theseus::Alignment alignment_object = aligner.align(sequence, start_vertex, start_offset, use_density_drop, use_lag_pruning);
```

If you only need the alignment score (and where the alignment ends), you can use the **score-only mode**. In this mode Theseus does not perform the backtrace, so it does not keep the information needed to recover the alignment, which reduces the memory usage. Enable it with:
```
aligner.set_alignment_mode(theseus::AlignmentMode::ScoreOnly);
```
The returned Alignment object contains the score (`alignment.score`), the start and end offsets, and a path with only the last vertex of the alignment. The CIGAR (`edit_op`) is empty. The score is the same as in the default mode (`theseus::AlignmentMode::Full`, used by default). Use `Full` again to recover the complete alignment. The score-only mode is only available for the sequence-to-graph aligner, the MSA aligner always computes the full alignment.

### <a name="graph_creation"></a> 2.3. Creating a graph

The Theseus' library, allows you to create your own reference graphs to perform sequence-to-graph alignment. A graph is composed of two key elements: nodes and edges. Nodes store genomics' data in the form of a sequence of characters, and edges represent connections between these existing nodes. If you want to create a graph, you first have to include the "theseus/graph.h" header file:
```
#include "theseus/graph.h"
```

Then, you have to assign it a name upon creation:
```
theseus::Graph graph;
```

You can add nodes using the "add_node(string)" functionality, that always returns a node identifier (NodeId) that will be necessary to create edges. For instance, if you want to create a node storing the sequence "ACTTAG" and another containing the sequence "T", you should:
```
theseus::Graph::NodeId id_node1 = graph.add_node("ACTTAG");
theseus::Graph::NodeId id_node2 = graph.add_node("T");
```

Now, if you want to create an edge connecting these two existing nodes, you have to use the "add_edge(id_node1, id_node2)" functionality:
```
graph.add_edge(id_node1, id_node2);
```


### <a name="penalties_usage"></a> 2.4. Penalties and dual-affine gaps

Theseus works with minimization penalties: the match penalty is the lowest one (0 by default, it can also be positive) and the rest are higher-valued. By default, gaps are scored with a single affine function, where a gap of length `k` costs `gap_open + k * gap_extend`:
```
theseus::Penalties penalties(match, mismatch, gap_open, gap_extend);
```

Theseus also supports **dual-affine gaps** (two-piece affine gaps). In this case the gap function has two pieces, and a gap of length `k` costs the cheapest of them:
```
min(gap_open + k * gap_extend, gap_open2 + k * gap_extend2)
```
The first piece is cheaper for short gaps (lower opening but higher extension), and the second one for long gaps (higher opening but lower extension). To use it, give the second piece to the `Penalties` constructor:
```
theseus::Penalties penalties(match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2);
```
The second piece must satisfy `gap_open2 >= gap_open` and `match < gap_extend2 <= gap_extend`, otherwise the aligner throws an `std::invalid_argument` exception. The penalties object is the only thing to change: the rest of the code (MSA, alignment to a graph, tools) is used as before.

In the command line tools (see [Tools](#theseus_tools)), the second piece is set with `-O/--gapo2` and `-E/--gape2`, which must be given together.

The score of an alignment according to these penalties is available in `alignment.score`.

## <a name="theseus_tools"></a> 3.Tools

The Theseus library implements two minimal tools to use the Theseus algorithm on the consensus generation context. It is important to note that these tools are not production ready.

### <a name="consensus_tool"></a> 3.1. MSA tool

This example illustrates how to use the **MSA** tool. This tool computes the MSA of the set of sequences in an given input *.fasta* file, allowing to add partial sequences, as long as they start on either end of a backbone sequence. The executable is located in the path */build/tools/theseus_msa*:
```
cd build/tools/
```

Select the scoring scheme, set the input and output files and execute the tool. Each execution of *theseus_msa* lets you select the following parameters:
```
Usage: theseus_msa [OPTIONS]
                 Options:\n"
                   -m, --match <int>           The match penalty                                       [default=0]
                   -x, --mismatch <int>        The mismatch penalty                                    [default=2]
                   -o, --gapo <int>            The gap open penalty                                    [default=3]
                   -e, --gape <int>            The gap extension penalty                               [default=1]
                   -O, --gapo2 <int>           The second gap open penalty (dual affine)               [default=none]
                   -E, --gape2 <int>           The second gap extension penalty (dual affine)          [default=none]
                   -t, --output_type <int>     The output format of the multiple alignment             [default=0=MSA]
                                                0: MSA: Standard Multiple Sequence Alignment format,
                                                1: GFA: Output the resulting POA graph in GFA format,
                                                2: Consensus - Heaviest Bundle,
                                                3: Consensus - Weighted Majority Voting,
                                                4: Dot: Output in .dot format. Only tractable for small graphs
                   -f, --output <file>         Output file                                             [Required]
                   -s, --sequences <file>      Dataset file                                            [Required]

                  Heuristics:
                   -l  --lag_pruning           Activate the pruning of diagonals lagging behind in the alignment.
```

An example of the execution of *theseus_msa* is shown in the following piece of code
```
./theseus_msa -m 0 -x 2 -o 3 -e 1 -t 0 -f output_file.out -s sequences.fasta
```


### <a name="seq_to_graph_tool"></a> 3.2. Seq-to-graph tool: theseus_aligner
This example illustrates how to use the **theseus_aligner** tool. This tool aligns a set of sequences, given their starting vertices, offsets, and orientations, to a reference graph. Two input files are required: 1) The reference graph in *.gfa* format, 2) The sequences to be aligned with the starting alignment positions and orientation in *.fasta* format.

**[IMPORTANT]**
The *.fasta* file containing sequences has a special structure. As all .fasta files, the data associated to each sequence has two parts: 1) A line starting with ">" containing metadata, and 2) the sequence itself, that appears on the next lines.

1) Constists of three elements: *start_vertex* name, *start_offset*, orientation of the alignment. The orientation of the alignment can either be "+" or "-" and it indicates Theseus whether to align from the forward or reverse strand of the given start node. *[IMPORTANT]*, the *start_offset* value is respect the starting node in the given alignment orientation. That is, if alignment orientation is reverse ("-") and the offset is 0, alignment starts at the beggining of the reverse complemented segment.
2) Contains the sequence itself.

For instance, a sequence ACGT starting at vertex *v1*, offset *3*, in the forward orientation *+*, and aligned from left to right would be codified as:
```
> v1 3 +
ACGT
> ... (following lines)
```

The executable is located in the path */build/tools/theseus_aligner*:
```
cd build/tools/
```

Select the scoring scheme, set the input and output files and execute the tool. Each execution of *pericles* lets you select the following parameters:
```
Usage: theseus_aligner [OPTIONS]
                 Options:
                   Penalties:
                   -m, --match <int>            The match penalty                                       [default=0]
                   -x, --mismatch <int>         The mismatch penalty                                    [default=2]
                   -o, --gapo <int>             The gap open penalty                                    [default=3]
                   -e, --gape <int>             The gap extension penalty                               [default=1]
                   -O, --gapo2 <int>            The second gap open penalty (dual affine)               [default=none]
                   -E, --gape2 <int>            The second gap extension penalty (dual affine)          [default=none]

                   I/O:
                   -g, --graph_file <file>      Graph file in .gfa format                               [Required]
                   -s, --sequences_file <file>  Sequences and starting positons in .fasta format        [Required]
                   -f, --output_file <file>     Output file                                             [Required]

                  Heuristics:
                   -d  --density_heuristic     Activate the drop heuristic based on advancement density.
                   -l  --lag_pruning           Activate the pruning of diagonals lagging behind int the alignment.
```

An example of the execution of *pericles* is shown in the following piece of code
```
./theseus_aligner -m 0 -x 2 -o 3 -e 1 -g reference_graph.gfa -s sequences.fasta -f output.out
```


## <a name="theseus_heuristics"></a> 4. HEURISTICS
Theseus library implements some heuristic approaches that accelerate alignment at the expense of a limited loss in accuracy. In particular, Theseus implements 1) a **pruning heuristic** that discards diagonals that have fallen behind in the alignment, as long as the alignment has shown a significant advancement in the last scores, and 2) a **drop heuristic** that drops alignment when the advancement density (number of offsets advanced in the last scores) is very low. You can activate these heuristics when calling the align functionality in **theseus::TheseusAligner** or a **theseus::TheseusMSA** aligners.

Moreover, the two minimal tools provided in this alignment library allow you to activate these two heuristics from the command line, so that they are applied to all alignments in the instantiated Aligner object. You can use both heuristics in the *theseus_aligner* tool, and the lag pruning heuristic in the *theseus_msa* tool.


## <a name="theseus_datasets"></a> 5. DATASETS

The datasets used in our paper are available in Zenodo. [Link to datasets](https://zenodo.org/records/18482097)


## <a name="theseus_bugs"></a> 6. REPORTING BUGS AND FEATURE REQUEST

Feedback and bug reporting is highly appreciated. Please report any issue or suggestion on github or email to the main developer (albert.jimenez1@bsc.es). Don't hesitate to contact us if:
  - You experience any bug or crash.
  - You want to request a feature or have any suggestion.
  - Your application using the library is running slower than it should or you expected.
  - Need help integrating the library into your tool.


## <a name="theseus_licence"></a> 7. LICENSE

Theseus-lib is distributed under MIT licence.


## <a name="theseus_authors"></a> 8. AUTHORS

Albert Jimenez Blanco (albert.jimenez.blanco@upc.edu) is the main developer and the person you should address your complaints.

Lorién López-Villellas has had major contributions both in the technical implementation of Theseus and the final structure of the library.

## <a name="theseus_cite"></a> 9. CITATION

If you use Theseus in your research, please cite:

**Albert Jiménez-Blanco, Lorién López-Villellas, Juan Carlos Moure, Miquel Moreto, Santiago Marco-Sola**. ["Theseus: fast and optimal affine-gap sequence-to-graph alignment"](https://academic.oup.com/bioinformatics/article/42/9/btag561/8742314). Bioinformatics, 42(9), btag561, 2026. DOI: [10.1093/bioinformatics/btag561](https://doi.org/10.1093/bioinformatics/btag561)

```bibtex
@article{jimenezblanco2026theseus,
  title   = {Theseus: fast and optimal affine-gap sequence-to-graph alignment},
  author  = {Jim{\'e}nez-Blanco, Albert and L{\'o}pez-Villellas, Lori{\'e}n and Moure, Juan Carlos and Moreto, Miquel and Marco-Sola, Santiago},
  journal = {Bioinformatics},
  volume  = {42},
  number  = {9},
  pages   = {btag561},
  year    = {2026},
  doi     = {10.1093/bioinformatics/btag561}
}
```

