### nfcore-mag-shiny

This web app is a graphical interface for the [nf-core/mag](https://nf-co.re/mag) pipeline. Check the nf-core/mag webpage for the meaning of the different options and their usage.

### Usage

 Select a folder containing fastq files using the `Select fastq folder` button and press `Run nf-core/mag pipeline`.

- By deafalt, results are written to the `results` folder within the fastq folder
- To test the pipeline - select nextflow profile `test` under more options (no need to select fastq folder in this case)
- You can always restart the app by clicking on the title
- You can start the run with Nextflow Tower, then go to the provided url and monitor it there!

The nextflow log files can be accessed by navigating to the fastq folder and running

```bash
cat .nextflow.log
```

### Questions

The Shiny app is written and maintained by [Angel Angelov](https://github.com/angelovangel), the nf-core/mag pipeline - by the [nf-core team](https://nf-co.re/mag).
