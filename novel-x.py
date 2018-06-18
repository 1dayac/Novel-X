import click
import pysam
import json
from os import mkdir, path, symlink, chdir
from shutil import copy2
from subprocess import call, Popen

def get_read_group(bam):
    readgroup = ""
    bamfile = pysam.AlignmentFile(bam, "rb")
    for read in bamfile.fetch(until_eof=True):
        try:
            readgroup = read.get_tag('RG')
            break
        except:
            pass
    bamfile.close()
    return readgroup.replace(":", "_")[:-2]

def create_config(bam, genome, nt, outdir, lr20):
    data = {}
    data['genome'] = genome
    data['additional_flags'] = "--lr20" if lr20 else ""
    data['blast_db'] = nt if nt != "" else "None"
    data['root'] = path.dirname(__file__)
    data['sample'] = path.basename(bam)[:-4]
    data['outdir'] = outdir
    data['readgroup'] = get_read_group(bam)
    with open(outdir + "/config.json", 'w') as configfile:
        json.dump(data, configfile, sort_keys=True, indent=4)



@click.group()
def main():
    pass

@main.command()
@click.option('--outdir', nargs = 1, required = True)
def restart(outdir):
    """Restart unfinished 10X-pipeline for novel insertion detection."""
    chdir(outdir)
    process = Popen(['snakemake', '--unlock'])
    process.wait()
    process = Popen(['snakemake'])
    process.wait()

@main.command()
@click.option('--bam', help = "No options defined but a name was passed", required = True)
@click.option('--genome', help = "Genome file in fasta or fasta.gz format", required = True)
@click.option('--nt', default = "", nargs = 1, help = 'Folder containing NT database. '
                                      'If not provided filtering of non-human sequences is not performed')
@click.option('--outdir', nargs = 1, required = True)
@click.option('--lr20', default=False, help = 'If your BAM-file was produced by LongRanger 2.0 you should provide '
                                              'this option to avoid failures')
def run(bam, genome, nt, outdir, lr20):
    """Run 10X-pipeline for novel insertion detection."""
    try:
        mkdir(outdir)
    except:
        print("Output folder can't be created")
        return -1
    create_config(bam, genome, nt, outdir, lr20)
    copy2(path.dirname(path.realpath(__file__)) + "/path_to_executables_config.json", outdir)
    copy2(path.dirname(path.realpath(__file__)) + "/Snakefile", outdir)
    mkdir(outdir + "/sample")
    symlink(bam, outdir + "/sample/" + path.basename(bam))
    chdir(outdir)
    process = Popen(['snakemake'])
    process.wait()

if __name__ == '__main__':
    main()