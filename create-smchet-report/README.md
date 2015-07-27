Description
===========
These scripts will generate an intermediate JSON representation from PhyloWGS'
tumor phylogenies. From these, it will then produce formats needed for the
[SMC-Het Challenge](http://dreamchallenges.org/project/home-upcoming/dream-9-5-icgc-tcga-dream-somatic-mutation-calling-tumor-heterogeneity-challenge-smc-het/).


Example usage
=============

First, write JSON files describing the dataset:

    PYTHONPATH=/path/to/phylowgs_dir \
      python2 write_tree_json.py \
      fd201154-d55b-4be0-a458-2eb92a9f2fe8 \
      /path/to/trees.zip \
      summ.json.gz \
      muts.json.gz \
      mutass.zip

Then, create the Challenge outputs based on the JSON representation:

    mkdir outputs
    python2 write_report.py \
      summ.json.gz \
      muts.json.gz \
      mutass.zip \
      outputs


Full usage
==========
JSON writer:

    usage: write_tree_json.py [-h] [--num-trees NUM_TREES]
                              dataset_name tree_file tree_summary_output
                              mutation_output mutation_assignment_output

    Write JSON files describing trees

    positional arguments:
      dataset_name          Name identifying dataset
      tree_file             File containing sampled trees
      tree_summary_output   Output file for JSON-formatted tree summaries
      mutation_output       Output file for JSON-formatted list of mutations
      mutation_assignment_output
                            Output file for JSON-formatted list of SSMs and CNVs
                            assigned to each subclone

    optional arguments:
      -h, --help            show this help message and exit
      --num-trees NUM_TREES, -n NUM_TREES
                            Only examine given number of trees (default: None)


Challenge-outputs writer:

    usage: write_report.py [-h] [--min-ssms MIN_SSMS]
                           tree_summary mutation_list mutation_assignment
                           output_dir

    Write SMC-Het Challenge outputs

    positional arguments:
      tree_summary         JSON-formatted tree summaries
      mutation_list        JSON-formatted list of mutations
      mutation_assignment  JSON-formatted list of SSMs and CNVs assigned to each
                           subclone
      output_dir           Directory in which to save Challenge outputs

    optional arguments:
      -h, --help           show this help message and exit
      --min-ssms MIN_SSMS  Minimum number of SSMs to retain a subclone (default: 3)

Outputs
=======
* `1A.txt`: Inferred cellularity of the sample
* `1B.txt`: Best guess for number of cancerous populations
* `1C.txt`: For each cluster, 3 columns: the cluster number (starting from 1),  the typical number of ssms assigned to it and the phi value of the cluster.  
* `2A.txt`: Assignment of SSMs to clusters. Row i indicates will list the cluster index for SSM number i (i.e., the SSM with the identifier `s<i - 1>`).
* `2B.txt`: Full NxN co-clustering matrix
* `3A.txt`: For each of the best guess clusters, 2 columns: cluster ID and the cluster ID of its parent (0 is root node)
* `3B.txt`: NxN Ancestor-decedent matrix. Entry i,j = The probability that i is in node that is an ancestor of node containing j. 
