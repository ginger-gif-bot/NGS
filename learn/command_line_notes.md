# Bash notes 
---
# Note
- `touch` makes new files 
- `mkdir` makes new folder

- `rm` deletes a file
- `rmdir` deletes a folder (cannot delete a folder if folder is not empty)
- `rm -r foldername` deletes a folder even if the folder is not empty


## Navigation
1. `pwd`= print working directory(where i am)

1. `ls`= shows all the files in that folder
    - `ls -lt` = shows extra info about the files
    - `ls -ltr` = Lists files sorted by time modified, showing the newest files last
    - `ls -lart` = shows the hidden files in a folder
    - `ls -lh` = Lists files with their details (permissions, date modified) and makes the file sizes human-readable

1. `cd <directory/foldername>`= change directory(move into folder)
    - `cd ..` = brings you back a step 
    - `cd /folderA/folder_inside_folderA` =  to go inside a folder's folder
    - `cd` / `cd /` = brings you to the root
    - `cd ../..` = brings two step back
    - `cd ~` = brings you to the home directory
    - `cd -` = brings to the previous folder you were in

## Useful Tricks
1. `whoami` = tells the user that is using it currently 

1. `date` / `date +%D(date)/T(time)/H(hour)`etc = to check the date and time

1. `clear` = to clear the terminal
    - `ctrl + l` = shortcut key to clear the terminal

1. `man bash_command` = gives a manual for that command

1. `chmod num_code file/foldername` = gives permisiion 
    - permission to owner, worker, public
    - the num_code can be calculated via chomod calculator online

1. `gzip filename` = to zip a file
    - `gunzip filename.gz` =  to unzip/extract a file
    - `zcat filename.gz` = to view a compressed file without unzipping it

## Creating, Copying, Moving, Deleting
1. `mkdir <folder/directoryname>` = make new folder

1. `rmdir <folder/directoryname>` = delete a folder

1. `cp path/filename_you_want_to_copy  path/filename_you_want_to_copy_in` = to copy directories
    - `cp -i path/filename_you_want_to_copy  path/filename_you_want_to_copy_in` = warns if you have the file already coz cp overwrites the files

1. `/mnt/c/Users/dell/Downloads` = to get a file from downloads

1. `vim filename` = creates/opens a file and lets you edit in the terminal itself
    - press `i` to start writing 
    - write your code
    - press `esc` then `:wq` to exit and save
    - press `esc` then `:q!` to discard the changes and save & exit 

1. `mv filename/path` = renames/moves your files
    - `mv previous_filename new_filename` = renames the file beacuse you are in the same folder
    - `mv new_filename folder/new_filename` = to move into new folder 
       (file should exist to move)

1. `more <filename>` = to read page by page 

1. `touch <filename>` = to make a file

## Navigating files
1. `cat <filename>` =  to see the contents of a file

1. `less <filename>` = to search contents in a file 
    - after doing this command a interacting like thing will open 
    - press `/`(forward slash) (to search from top to bottom) 
    - press `?` (question) (to search from bottom to top) 
    - type what you want to search 
    - press `n` to see other values with same spelling 
    - when you do less command, you are at the top of the file 
    - to get to the end of the file press `shift + g` 
    - to get to the top of the file press `p`

1. `head filename` = to see top 10 lines
    - `head -n 20` = to see top 20 lines, you can change the number

1. `wc` = finds lines, words, character counts
    - `wc -l` = counts lines
    - `wc -w` = counts words
    - `wc -c` = counts characters

1. `tail filename` = to see bottom 10 lines
    - `tail -n 20` = to see bottom 20 lines, you can change the number

## Searching and Filtering
1. `grep "word" filename` = to find word in a file (case sensitive)
    - `grep -i "word" filename` = to find word in a file (case *insensitive*)
    - `grep -c "word" filename` =  counts how many lines contained the word
    - `grep -v "word" filename` =  shows lines that do not contain the word
    * `grep ">" filename.fasta` = finds all the headers in fasta file
    - `grep -n ">"/"word" filename.fasta` = finds all the headers/words and also tells their line number in fasta file

1. `cat file_you_want_to_copy_from >> file_you_want_to_copy_in` = to copy the contents of one file to another without overwritting
    - `cat file1 file2 > combined_file.fasta` =  to combine contents of files into a third
    - `head/tail -n num file_you_want_to_copy_from > file_you_want_to_copy_in` =  to copy specific lines from onefile to other

1. `cut -f1 filename.tsv` = to extract a column from a tab separated file(you can change the num to see different columns)

1. `cut -d -f1 filename.csv` = to extract a column from a comma separated file(you can change the num to see different columns)

1. `nl filename` = to see the numberlines of the dataset

1. `sort filename` = sorts the lines in alphabetical order
    - `sort -n filename`= sorts numerically (acsending order)
    - `sort-rn filename` = sorts numerically (descsending order)

1. `uniq filename` = removes duplicates
    - `uniq -c` = counts how many times each line appears

## Pipe
- used to chain the commands
- the result of left command is give to the righ command as input and chanining continues
- passes data from command to command

## Redirects
- passes data between a command and a file.

1. `The Overwrite Redirect (>)`
    - This takes the output of a command and dumps it into a file. If the file already exists, it wipes out everything inside that file first and replaces it with the new data. If the file doesn't exist, it creates a new one.

1. `The Append Redirect (>>)`
    - This is the safer cousin of >. Instead of wiping out the destination file, it gently glues the new output onto the very bottom of the existing file.

1. `The Input Redirect (<)`
    - It feeds the contents of the file as input to the command, completely bypassing the keyboard. The command processes that data exactly as if you had manually typed it out line by line, character by character.

1. `The Error Redirect (2>)` 
    - It shows the error

## Installing Tools and Downloading dataset 
1. `sudo apt update` = Refresh the list of available Linux packages
1. `sudo apt install toolname` =  Install a Linux tool (e.g. samtools, fastqc)
1. `conda install toolname` = Install a bioinformatics tool via Conda
1. `conda install -c bioconda toolname` = Install from the bioconda channel (most bio tools)
1. `pip install packagename ` = Install a Python library (e.g. pandas, biopython)
1. `wget URL` =  Download a file from the internet
1. `wget -O myfile.gz URL` =  Download and save with a specific filename
1. `curl -O URL` =  Alternative to wget for downloading files

## Keyboard Shortcuts and Useful Tricks
1. `Tab` =  Auto-complete a command or filename — use this constantly!
1. `Up/Down arrow keys` =  Scroll through previous commands
1. `Ctrl + C` = STOP a running command immediately
1. `Ctrl + L` = or clear Clear the terminal screen
1. `Ctrl + A` = Jump to beginning of current line
1. `Ctrl + E` = Jump to end of current line
1. `man commandname` = Open the manual/help for any command (press q to exit)
1. `commandname --help` = Quick help for a command
1. `history` = Show list of all previous commands you have typed
1. `!!` = Repeat the last command

## awk
- like the if in python
- `awk 'if condition {then command}' filename` 
- `awk '$3 >30 {print $1 $2}' filename` = means if values in column 3 are above 30 then print column 1 and 2 for that
- `$0` = the entire line
- `$1` = column 1 and so on
- `NF` = total number of columns
- `NR` = current row records

### basename
- it strips away the directort path off a file, leaving just the filename
- `variablename=(basename "$iterator/var_containing_the_file_path/name")` 
- use the variablename after that instead of iterator name

### To initialize  a bash file
- `#!/bin/bash`
- at the top of the script
- to run= `bash scriptfile`

###
`$(...)` -> command substitution
`$variable` -> variable expansion
`if [ -f... ]` -> to check if file exists
##### To replace 
variable=`"$(var_you_want_to_perform_on/thing_to_replace/thing_to_replace_into)"`
##### other replacing way:
varr=`$(variable%.fastq)` -> means remove fastq from the end of the filename