# What to do if I want to launch a code on cluster?

No worries, my friend, here is a simple instruction! 

First, connect to the cluster by typing the following command in your terminal:
```shell
[yggdrasil@sowelu ~]$ ssh gene@login.mesocentre.univ-amu.fr
```

Type the password if prompted. Done, now you are inside a shell! Now you will have to change your working directory to the one where all our files are stored. 

```shell
[gene@login01 ~]$ cd /scratch/gene/
[gene@login01 gene]$ ls
```

Now you can see all the folders containing in our working directory! Cool, no? 

## What, we have to use the terminal? 

Terminals are fun! Here is a small reminder for the usage of the termminal which you *may* need during this session. For the longer list you should better google something like "Terminal Cheat Sheet". So,

* **pwd** - prints your working directory
* **ls** - lists all the files in the directory. For a more detailed list of options use **ls -il**
* **cd** - changes the working directory
* **rm** - removes the unnecessary files
* **cat** - shows the content of the text file in terminal. If the file is long, it may be better to use **tail** or **head** instead. 

Now, the script, which launch the cluster, is called *slurmjob.sh*. But before launching the script, you have to put your own script to cluster in the folder *horizontal_gene_transfer/Mesocentre/*. It is very simple! Put the code from your machine to your GitHub folder with the same name, and then push your changes to the server. It is done! Now, go inside *horizontal_gene_transfer* (using **cd**, of course) and pull the changes. 

## SLURM TIME!

Now, when everything is settled, it is time for a game! Before redacting the slurm script, please check if your Python code includes *matplotlib.use("Agg")*, as otherwise it will throw an error when calling a matplolib! Also check if includes the right path to the directories, so that all the figures are stored in */scratch/gene/horizontal_gene_transfer/Mesocentre/Figures*. 

Finally, change the last command in *slurmjob.sh*. It must include a path to your script. You have to do that with vim (yes-yes!). Check the short tutorial in links to this page. Or, if you feel like it, you can write your own script inside a GitHub folder, copypasting the content of the one in */scratch/gene/*, or the one from the tutorial on Aix-Marseille University website (which is equivalent), and push your changes to the folder, and launch it. 

You can do it and see how it works with the following commands:

```shell
[gene@login01 gene]$ sbatch ./slurmjob.sh
[gene@login01 gene]$ squeue -u gene
[gene@login01 gene]$ ls -il logs/
```

If everything is fine, you will see your job ID in queue (second command), and two new files in folder *logs* with extension *.out* and *.err*. The last must be empty, if everything is fine. If you don't see any jobs after the second command - check the error file. The job has either finished, or failed to launch because, for example, a mistake in your code. 

You can also observe how your job is going (if, for example, you print a message every 100-th iteration) by printing the changes in the log file in your console, for that you may use:
```shell
[gene@login01 gene]$ tail -f logs/your_log_file.out
```
**Tips&Tricks:** If you think that your code is not writing anything, while it should, be sure to add an additional argument **flush=True** in your Python's **print** command. 

To stop the job and to disconnect from mesocentre (you are not required to cancel the job upon disconnecting, of course):

```shell
[gene@login01 gene]$ scancel ID_OF_THE_JOB
[gene@login01 gene]$ exit
```

To see the preliminary results (for example, some messages which you usually print terminal, check the folder *logs*, it contains the files with an extension *.out*, where you can see the text output of your code. 

Now, when you got what you want (for example, new figures, or some new data, or whatever), push your changes to the server using *git* command inside a *horizontal_gene_transfer* folder (as usual). 

If something does not work or you are just scared by a colorscheme in a terminal - contact me (Anya), I will do my best to help you! If I am not available, and/or if the problem is too sophisticated - contact Benoit. If you want to fix everything by your own means, use google or the other links from this list:
* [How to work in Vim](https://vim.rtorr.com)
* [How to use Marseille's mesocentre](https://mesocentre.univ-amu.fr/tutoriauxn/)


Have fun! 
