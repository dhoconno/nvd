<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>

<!-- ABOUT THE PROJECT -->
## About NVD
![NVD flowchart](nvd_flow.png)

NVD is an opinionated workflow for identifying and exploring metagenomes for families of viruses that infect humans. It is opinionated in that it focuses narrowly on these virus families to the exclusion of other potentially interesting viral and non-viral microbial taxa. NVD works with arbitrarily large Oxford Nanopore and Illumina datasets without downsampling. The results are available in a LabKey Server data explorer available to O'Connor lab users and collaborators.

### Why Another Metagenomic Workflow?

Good question. When I started working on NVD in June 2024, I hoped to find a widely used, easy-to-use workflow optimized for viral metagenomics. CZID is excellent, but downsamples to a subset of reads. Sometimes this isn't a problem; sometimes it causes substantial reductions in sensitivity. The outstanding Sourmash tool performs informative, efficient searches against reference databases but doesn't generate contigs that can be used to characterize interesting viruses.

### What Does NVD Do?

NVD begins with an Illumina paired-end or Oxford Nanopore dataset. It can also evaluate datasets on NCBI SRA. The file integrity of these reads is checked and, in the case of Illumina datasets, paired-end sequences are interleaved and properly paired using the bbmap `repair.sh` tool.

These reads are the classified with the NCBI STAT tool, specifically the `aligns_to` command. A tax list containing the entire subtree of the virus families that infect humans is used for classification. This means that NVD is designed to find members of these virus families:

- Adenoviridae
- Anelloviridae
- Arenaviridae
- Arteriviridae 
- Astroviridae
- Bornaviridae
- Peribunyaviridae
- Caliciviridae
- Coronaviridae
- Filoviridae
- Flaviviridae
- Hepadnaviridae
- Hepeviridae
- Orthoherpesviridae
- Orthomyxoviridae
- Papillomaviridae
- Paramyxoviridae
- Parvoviridae
- Picobirnaviridae
- Picornaviridae
- Pneumoviridae
- Polyomaviridae
- Poxviridae
- Sedoreoviridae
- Retroviridae
- Rhabdoviridae
- Togaviridae
- Kolmioviridae

Reads that have kmers that match one or more viruses within these families are saved by `aligns_to.` These are "putative viral reads." NVD uses seqkit grep to extract these sequences from the original dataset.

The putative viral reads are then assembled into contigs with SPAdes. Empirically, the `--sewage` assembly preset provides the best results with both ONT and Illumina data, so this mode is used by NVD.

Contigs produced by SPAdes are then quality checked. Contigs with unusually high entropy, often consisting entirely of repetitive sequence, are filtered using `bbmask.sh.` Contigs shorter than 200bp after masking are removed, since short contigs are minimally informative. The filtered contigs are then classified with the NCBI STAT tool against a "coarse" reference library that is broadly representative of all the sequences in NCBI. A subtree of possible exact classifications is derived from the coarse classifications and used for more accurate classification. From these accurate classifications, contigs that appear to be viral are extracted.

These viral contigs are then used as BLAST queries against the NCBI `core-nt` database. Megablast classification typically identifies matches for most contigs. Those contigs that do not have megablast matches are then used as queries for blastn. 

Additionally, the putative viral reads are mapped against the viral contigs. This allows further exploration and confirmation of any unexpected or unusual hits.

The blast classifications, read mapping information, and associated contig metadata are then saved to an output folder in Excel format (and loaded into a database in the O'Connor lab LabKey Server, if appropriate access credentials are provided. `.zst` archives containing useful data files (e.g., read mapping `.bam` files) are also saved in the output folder.
  
<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- GETTING STARTED -->
## Getting Started

NVD is written in Snakemake and has its dependencies bundled in an Apptainer container. NVD should work on machines with at least 32GB of available RAM. The `resources.zst` file is approximately 230GB. After decompression, it is approximately 300GB. The amount of storage needed depends on the size of the FASTQ files; 530GB plus 3x the size of the input FASTQ files is a good guide. For example, if the input FASTQ files are 100GB, ensure that at least 830GB of storage are available.

### Prerequisites

- [Miniconda](https://docs.anaconda.com/miniconda/miniconda-install/)
- [Apptainer](https://apptainer.org/docs/admin/main/installation.html)
- [Snakemake](https://snakemake.readthedocs.io/en/stable/)
- The `resources.zst` archive containing NCBI BLAST `core-nt`, NCBI STAT databases, a taxonomic rank database, and a taxonomic list of the subtree of human-infecting virus families. Get this from DHO.
- The `workflow` folder containing the Snakemake workflow and associated python scripts
- A `config` folder containing a `config.yaml` file specifying runtime variables and the samples to be analyzed. You also need to get the LabKey API Key, LabKey username, and LabKey password information from DHO for LabKey integration.

### Installation

1. Create a conda environment with Apptainer and Snakemake (Linux)
	```
	conda create -n nvd
	conda activate nvd
	conda install -y snakemake apptainer -c conda-forge
	```
2. [Download](https://g-2e5b4e.dtn.globus.wisc.edu/nvd/nvd_30572.sif) the Apptainer image file.
3. Clone the repo to the working directory
   ```sh
   git clone https://github.com/dhoconno/nvd.git
   ```
4. [Download](https://g-2e5b4e.dtn.globus.wisc.edu/nvd/resources.zst) the `resources.zst` file containing databases and taxonomy files
5. Start an Apptainer shell in the working directory with the repo files
   ```
   apptainer shell nvd_30572.sif
   ```
6. Decompress the `resources.zst` file in the working directory
   ```
   tar -I zstd -xvf resources.zst
   ```
7. Copy gzipped-FASTQ files to process into `data` folder within the working directory.
8. Modify the `config.yaml` file to specify the samples to process and the path(s) to their FASTQ files. Here are example entries for the three supported file types:
   ```
   - name: water_S80_
    r1_fastq: data/water_S80_L006_R1_001.fastq.gz
    r2_fastq: data/water_S80_L006_R2_001.fastq.gz
   - name: AE0000100C7A40
    sra: SRR24010780
   - name: AE0000100A8B3C
    ont: data/AE0000100A8B3C.fastq.gz
   ```

After installation, there should be `data`, `config`, `resources`, and `workflow` folders in the working directory.
<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- USAGE EXAMPLES -->
## Usage

Use this space to show useful examples of how a project can be used. Additional screenshots, code examples and demos work well in this space. You may also link to more resources.

_For more examples, please refer to the [Documentation](https://example.com)_

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- ROADMAP -->
## Roadmap

- [ ] Feature 1
- [ ] Feature 2
- [ ] Feature 3
    - [ ] Nested Feature

See the [open issues](https://github.com/github_username/repo_name/issues) for a full list of proposed features (and known issues).

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to learn, inspire, and create. Any contributions you make are **greatly appreciated**.

If you have a suggestion that would make this better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement".
Don't forget to give the project a star! Thanks again!

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

<p align="right">(<a href="#readme-top">back to top</a>)</p>

### Top contributors:

<a href="https://github.com/github_username/repo_name/graphs/contributors">
  <img src="https://contrib.rocks/image?repo=github_username/repo_name" alt="contrib.rocks image" />
</a>



<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- CONTACT -->
## Contact

Your Name - [@twitter_handle](https://twitter.com/twitter_handle) - email@email_client.com

Project Link: [https://github.com/github_username/repo_name](https://github.com/github_username/repo_name)

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- ACKNOWLEDGMENTS -->
## Acknowledgments

* []()
* []()
* []()

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/github_username/repo_name.svg?style=for-the-badge
