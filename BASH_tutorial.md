# BASH introduction
In this repository, we provide toy data to learn the basic principles of bash.
In the toy folder, you can find files and subfolders to practice on.
This tutorial was initially developed for the GCRF STAR course, and can be found [here](https://github.com/evotools/GCRF_tutorial).

## Download the data
First, we download the toy data for the examples below using `wget`:
```
wget https://www.dropbox.com/s/7zkjy2nt55fdopk/toy.zip
```

Then, we unzip the folder using the `unzip` command:
```
unzip toy.zip
```

...and then enter the unzipped folder:
```
cd toy/
```
From here on, we can start learning how to move around in bash!

## Directories
To see the current working directory, use the command:
```
pwd
```

To list the content of the current directory, use:
```
ls
```

If we want more details about the content, we can use the following options:
```
ls -alh
```
This will show a detailed list of files and folders, including information on:
1. Their type (d) and permissions
2. The owner and the group
3. The file size 
4. When it has been created
5. The folder name 

If we want to see the content of a different folder, I can simply specify its path:
```
ls /home
```

## Changing directory
To change directory, I can simply use the `cd` command:
```
cd data
pwd
```
To move to the parent directory, I can use the `../` path:
```
cd ../
pwd
```

## Creating and removing an empty directory
To create a directory, I can use the `mkdir` command:
```
mkdir newdir
```
And we can see the new folder with `ls`:
```
ls
```

I can then remove it using the `rmdir` command:
```
rmdir newdir
```
## Creating and deleting an empty file
To create a directory, I can use the `touch` command:
```
touch myfile.txt
```

If I use the command on an existing file, it will cause it to update the creation date
to the moment the command is run:
```
ls -alh myfile.txt
```
Shows `-rwxrwxrwx 1 user user    0 May 23 13:00 myfile.txt`
If I touch it again, it will update the date:
```
touch myfile.txt
ls -alh
```
Will show `-rwxrwxrwx 1 atalent atalent    0 May 17 13:01 myfile.txt`.

I can then remove the file using the `rm` command:
```
rm myfile.txt
```

## Copy, move and rename files
To create a copy of a file, simply use the `cp` command:
```
cp chromosomes.txt cow_chromosomes.txt
```
As you can see, my copy has a different name than the original. 

I can copy to another location simply by specifying the full or relative destination path:
```
cp chromosomes.txt ./data/
```
This will create a copy of `chromosomes.txt` in the `data` subfolder. 
And see the new file in `data` with `ls` again:
```
ls ./data
```

I can also rename while copying the file simply by specifying the new name in the command as a 
second argument:
```
cp chromosomes.txt ./data/cow_chromosomes.txt
```

If I want to move a file from one directory to another, I can use the `mv` command.
To try this command out, let's first create a new folder:
```
mkdir files
```
Then, we can move the `cow_chromosomes.txt` from the `data` folder to our second folder:
```
mv ./data/cow_chromosomes.txt ./files
```
The syntax is quite simple: move `cow_chromosomes.txt` from the `data` folder to the `doc` folder.
We can also rename the file simply by specifying the new name:
```
mv ./files/cow_chromosomes.txt ./files/cattle_chromosomes.txt
```
This will rename the file from `cow_chromosomes.txt` to `cattle_chromosomes.txt`.

## See the content of a file
We can see the content using a series of possible commands. 

We can see the first few lines of a file using `head`:
```
head ./data/toy.log
```
or the last few lines using `tail`:
```
tail ./data/toy.log
```
It is possible to specify the number of lines with the `-` followed by the number of lines:
```
head -15 ./data/toy.log
```
This will show the first 15 lines of the log file.

We can print the content of a file using `cat`:
```
cat ./data/toy.log
```
If we redirect (`>`) the output to a new file, we can save a copy:
```
cat ./data/toy.log > ./files/toy.doc2
```
If we provide two files, we can concatenate them:
```
cat ./data/toy.log ./files/toy.doc2 > ./files/toy.doc3
cat ./files/toy.doc3
```
If a file is compressed with gzip (ending with `.gz`), we can use zcat:
```
zcat ./data/toy.vcf.gz
```

Alternatively, we can see the content one screen at a time using `more`:
```
more ./data/toy.log
```
With this command, you can move to the next screen with the `enter` key, and close with the `q` key.
An alternative to `more` is `less`:
```
less ./data/toy.vcf.gz
```
Unlike `more`, `less` lets you move backwards in the file using the up/down arrows.
It also allows users to view the content of gzip-compressed files, such as our vcf file.
To close, simply press the `q` key.

https://github.com/DomeJoyce/epicrops-course-2026
## Looking for terms/words in a file
We can screen for words in a file using `grep` (**g**eneral **r**egular **e**xpression **p**rint):
```
grep variant ./data/toy.log
```
This will search for all partial occurrences of variant (i.e. "variant", "variants", "invariant").
If we want to seek *exact* matches, we can add the `-w` option:
```
grep -w variants ./data/toy.log
```

## A look in the past
We might need to re-run commands we recently executed. To facilitate this operation and avoid
rewriting everything, bash lets us see the past commands by using the `history` tool:
```
history
```
This will print a list of all commands we executed, anticipated by a unique number for the command:
```
history
 2009  ls
 2010  pwd
 2011  ls -al
 2012  history
```

At this point, we can proceed in two ways: 
1. Copy and paste the command: in putty, highlight the command, and then right-click twice
2. Refer to the number 

The second can be more convenient, since it saves us from copying and pasting. We can simply refer to a command in the history
by using the `!` symbol, followed by the number as follows:
```
!2011
```
In this case, using `!2011` is the same as typing `ls -al`.
Remember that whenever bash closes, the history is refreshed, and the number won't necessarily match.

# Special characters
Characters with different behaviour in the BASH console are considered "special".
These characters allow the user to control the flow of the operations and combine multiple of them.

The first special character is the new line character (coded as `\n`). This character signals to bash
that the command input is completed, and to try to run it.

Another important special character is the space (` `), which indicates the separation between commands 
or parameters. For this reason, it is better to avoid using it when naming a file, or replace it with 
other characters (e.g. `.`, `-` or `_` generally works well).

We have also already seen the `!` symbol, which can be used to recall specific commands from the history.

Other major special characters we will discuss here with examples are:
1. `*`
2. `&`
3. `\` (escape)
4. `|` (pipe)
5. `#` (comment)
6. `>` (redirection)
A more comprehensive list of characters can be seen well documented [online](https://www.cs.ait.ac.th/~on/O/oreilly/unix/upt/ch08_19.htm) 

## The * character
The first important character to know is the `*` character, which means "all". For example, using the command:
```
ls ./data/*
```
Will list every file in a folder. Changing it to:
```
ls ./data/*.log
```
Will list all files ending with the `.log` suffix. Alternatively, using the command:
```
ls ./data/toy*
```
Will list all files starting with the `toy` name. 

## The & character
The second special character is the `&` character. This character can have two separate meanings:
1. When it is used singularly at the end of a command, it runs the command in the background
```
cp ./data/toy.vcf.gz ./data/toy.2.vcf.gz &
```
2. When two are placed at the end of a command, they allow to run a second command IF the first is successful.
```
cp ./data/toy.vcf.gz ./data/toy.2.vcf.gz && mv ./data/toy.2.vcf.gz ./data/mytoy.vcf.gz
```
This command will try to copy `toy.vcf.gz` and, if this is successful, it will rename it.

## The \ (escape) character
The `\` character is called the 'escape' character and allows you to ignore the special effect of a character.
For instance, using the command:
```
cp ./data/toy.vcf.gz ./data/toy\ data.vcf.gz 
```
This will copy our `toy.vcf.gz` and rename it as `toy data.vcf.gz`: the newly created file will contain a space 
in the name. If I want to refer to the file in my next commands, I have to add `\` just before the space, in order
to allow the bash to correctly identify it.

Another common use is, for example, to break up a complex command into multiple, simpler lines:
```
cp ./data/toy.vcf.gz ./data/toy.2.vcf.gz && \
    mv ./data/toy.2.vcf.gz ./data/mytoy.vcf.gz
```
In the above example, the `\` after the `&&` allows us to consider the new line as a simple character and that the command is not
complete.

## The | (pipe) character
The `|` character allows moving the output of a command into the next one.
For example, if we want to view line 10 of a compressed file, I can combine the `zcat` command with the 
`head` command:
```
zcat ./data/toy.vcf.gz | head
```
Another useful example is how to view a specific line in a file. I can do this, for example, by combining the `head` and `tail` 
commands. For example, I can extract only line number 13 from my compressed vcf file as follow:
```
zcat ./data/toy.vcf.gz | head -13 | tail -1
```
In order, this command will:
1. Print to screen the content of the file
2. Extract the first 13 lines from the file
3. Show only the last one of them (line 13)

## The # (comment) character
The `#` character signals to the console that the text following it is a **comment**. Comments are lines, or parts of lines, 
that are ignored by the console.
For instance, the commands
```
# zcat ./data/toy.vcf.gz | head -13 | tail -1
```
Will cause nothing to happen. 

On the other hand, using the command:
```
zcat ./data/toy.vcf.gz | head -13 #| tail -1
```
Will cause the first two commands (`zcat` and `head`) to work, but not the last (`tail`). 

## The > (redirection) character
The `>` character signals to the console that the output from a command is sent from the screen (known as `stdout`) to a 
text file with a given name.
For instance, this command prints "Hello everyone!" to the screen:
```
echo "Hello everyone!"
```
But if we add the `>` followed by a file name **after** the command, it will save to a file:
```
echo "Hello everyone!" > Hello.txt
```
We can see the content of `Hello.txt` using the `cat` command:
```
cat Hello.txt
```
If a file already exists, the command will **overwrite** the content (use it with care!):
```
echo "Hello guys!" > Hello.txt
cat Hello.txt
```

If we use two `>` symbols (`>>`), it is possible to add lines to a file without overwriting:
```
echo "Bye guys!" >> Hello.txt
```
Then, we can see the content of the file with `cat` once more:
```
cat Hello.txt
```
