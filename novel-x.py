import click
from click_option_group import optgroup, MutuallyExclusiveOptionGroup
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

def create_config(bam, genome, nt, outdir, lr20, m, t, highcoverage, lowcoverage, tenx):
    data = {}
    data['genome'] = path.abspath(genome)
    data['additional_flags'] = "--lr20" if lr20 else ""
    data['blast_db'] = nt + "/nt" if nt != "" else "None"
    data['root'] = path.dirname(path.realpath(__file__))
    data['sample'] = path.basename(bam)[:-4]
    data['outdir'] = outdir
    data['threads'] = int(t)
    data['memory'] = int(m)
    data['memory_per_thread'] = int(m/(2*t)) + 1
    data['velvet_coverage'] = 8 if highcoverage else 2
    data['velvet_k_assembly'] = 63 if tenx else 49
    data['spades_k_assembly'] = 77 if tenx else 49
    data['tenx'] = "10X" if tenx else "other"
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
    process = Popen(['snakemake', '--unlock', '--cores', 'all'])
    process.wait()
    process = Popen(['snakemake', '--cores', 'all'])
    process.wait()

@main.command()
@click.option('--bam', help = "No options defined but a name was passed", required = True)
@click.option('--genome', help = "Genome file in fasta or fasta.gz format", required = True)
@click.option('--nt', default = "", nargs = 1, help = 'Folder containing NT database. '
                                      'If not provided filtering of non-human sequences is not performed')
@click.option('--outdir', nargs = 1, required = True)
@click.option('--lr20', default=False, is_flag = True, help = 'If your BAM-file was produced by LongRanger 2.0 you should provide '
                                              'this option to avoid failures')
@click.option('-m', default = 100, nargs = 1, help = 'Available memory specified in gygabytes')
@click.option('-t', default = 8, nargs = 1, help = 'Number of threads')
@optgroup.group('Coverage option group', cls=MutuallyExclusiveOptionGroup)
@optgroup.option('--high-coverage', 'highcoverage',  is_flag=True, default=True, help = "Flag for high-coverage data (60X or higher, default)")
@optgroup.option('--low-coverage', 'lowcoverage', is_flag = True, default = False, help = "Flag for low-coverage data (20-40X)")
@optgroup.group('Data type option group', cls=MutuallyExclusiveOptionGroup)
@optgroup.option('--10x', 'tenx', is_flag=True, default=True, help = "Default")
@optgroup.option('--stlfr', is_flag = True, default = False)
@optgroup.option('--tellseq', is_flag = True, default = False)
def run(bam, genome, nt, outdir, lr20, m, t, highcoverage, lowcoverage, tenx, stlfr, tellseq):
    """Run 10X-pipeline for novel insertion detection."""
    try:
        mkdir(outdir)
    except:
        print("Output folder can't be created. Probably it already exists")
        return -1
    if stlfr or tellseq:
        tenx = False
    if lowcoverage:
        highcoverage = False
    create_config(bam, genome, nt, outdir, lr20, m, t, highcoverage, lowcoverage, tenx)
    copy2(path.dirname(path.realpath(__file__)) + "/path_to_executables_config.json", outdir)
    copy2(path.dirname(path.realpath(__file__)) + "/Snakefile", outdir)
    mkdir(outdir + "/sample")
    symlink(path.abspath(bam), outdir + "/sample/" + path.basename(bam))
    chdir(outdir)
    process = Popen(['snakemake', '--cores', str(t)])
    process.wait()

if __name__ == '__main__':
    main()
