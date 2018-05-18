import click
import pysam

def get_read_group():
    pass

def get_base_folder():
    pass

@click.group()
def main():
    pass

@main.command()
@click.option('--outdir', nargs = 1, required = True)
def restart(outdir):
    """Restart unfinished 10X-pipeline for novel insertion detection."""
    pass

@main.command()
@click.option('--bam', type=click.File('rb'), help = "No options defined but a name was passed", required = True)
@click.option('--genome', type=click.File('rb'), help = "Genome file in fasta or fasta.gz format", required = True)
@click.option('--nt', nargs = 1, help = 'Folder containing NT database. '
                                      'If not provided filtering of non-human sequences is not performed')
@click.option('--outdir', nargs = 1, required = True)
@click.option('--lr20', default=False, help = 'If your BAM-file was produced by LongRanger 2.0 you should provide '
                                              'this option to avoid failures')
def run(bam, genome, nt, outdir, lr20):
    """Run 10X-pipeline for novel insertion detection."""
    pass

if __name__ == '__main__':
    main()