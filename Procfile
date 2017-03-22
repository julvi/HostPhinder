#################### How to run HostPhinder command line ######################
# NB: This version of HostPhinder has to be run From folder 
# /home/projects/pr_phage/people/juliav/HostPhinder_general/git_HP_repo/HostPhinder

./scripts/wrap_hostPhinder_cl.sh -i <input phage genome> -o <output directory> \
-t: <taxonomy [species, genus]>

#################### Run HostPhinder on queueing system ######################
xmsub -W group_list=pr_phage -A pr_phage \
-l nodes=1:ppn=2,mem=10gb,walltime=10:00 \
-d /home/projects/pr_phage/people/juliav/HostPhinder_general/git_HP_repo/HostPhinder \
-N <job_name> \ # Choose job name
-de \
-e <path/to/error/filename/errorfilename> \ # Choose error filename + path
-o <path/to/output/filename/outputfilename> \ # Choose output filename + path
-r y scripts/wrap_hostPhinder_cl.sh \
-i <input_fasta>\ # Input fasta file + path
-t species \ # Choose the taxonomy level of the prediction
-o <outputdir> # change to desired output folder + path

# Example on how to run as juliav
xmsub -W group_list=pr_phage -A pr_phage \
-l nodes=1:ppn=2,mem=10gb,walltime=10:00 \
-d /home/projects/pr_phage/people/juliav/HostPhinder_general/git_HP_repo/HostPhinder \
-N D9 \
-de \
-e /home/projects/pr_phage/people/juliav/HostPhinder_general/git_HP_repo/HostPhinder/D9.err \
-o /home/projects/pr_phage/people/juliav/HostPhinder_general/git_HP_repo/HostPhinder/D9.out \
-r y scripts/wrap_hostPhinder_cl.sh \
-i /home/projects/pr_phage/people/mvla/georgian_phages/data/D9_allContigs.fasta \
-t species \
-o /home/projects/pr_phage/people/juliav/HostPhinder_general/git_HP_repo/HostPhinder 

