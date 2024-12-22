* [Day1: Introduction to Unix/Linux OS](#day1-introduction-to-unix/linux-os)

  - [Introduction to Linux/Unix OS and command line](#introduction-to-linuxunix-os-and-command-line)
  - [Manipulating files, useful commands and tips](#manipulating-files-useful-commands-and-tips)
  - [Permissions, groups and control](#permissions-groups-and-control)
  - [Application of Linux](#application-of-linux)
  - [Environment variable](#environment-variable)
  - [Shell Scripting](#shell-scripting)
  - [Controlling tasks](#controlling-tasks)
  - [SSH into remote machine](#ssh-into-remote-machine)


# Day1: Introduction to Unix/Linux OS
## Introduction to Linux/Unix OS and command line
### What is Linux?
- UNIX is an Operating System (OS) initially developed in the 1960s.
- There are many different versions of UNIX, that share common similarities.
- The most popular varieties of UNIX are Solaris, Linux, and MacOS.
- UNIX systems have a graphical user interface (GUI) making them easier to use.
### Linux vs Unix
- Linux is a “clone” of the original Unix but doesn’t contain its code.
- Linux is free and open source, the original Unix is not (but some of its derivatives are).
- All command lines work the same on both.
### Why Linux?
- Linux is free and the most popular distributions are
Ubuntu, Fedora/Red Hat, Mandriva, etc.
- Low cost and very stable system.
- Most secure OS.
- Best multi-user and multi tasking OS.
- The world’s fastest super computers run Linux.
- Fast developing OS (many developers)
### Linux Distributions
Different Linux distributions are available. You can explore them here: [DistroWatch](http://distrowatch.com/).

- **Ubuntu**: A distribution that is easy and convenient to use for beginners.
- A simple guide to install Ubuntu on your machine can be found here: [Install Ubuntu Desktop](http://www.ubuntu.com/download/desktop/install-ubuntudesktop).

<p align="center">
<img src="https://github.com/Mo7ammedFarahat/MASRI-Jan25/blob/main/images/unix.png?raw=true" width="350"/>
</p>

<p align="center">
<img src="https://github.com/human-pangenomics/hprc-tutorials/blob/GA-workshop/assembly/genomics_aotearoa/images/sequencing/meryl_venn.png?raw=true" width="350"/>
</p>

### The Terminal
- A terminal refers to a wrapper program which runs a shell.
- There are many different Unix shells, the most popular shell for interactive use includes Bash: the default on most Linux installations.

<p align="center">
<img src="https://github.com/Mo7ammedFarahat/MASRI-Jan25/blob/main/images/Terminal.png?raw=true" width="650"/>
</p>
Even though it is a command line interface, the mouse is still handy (scroll, copy, paste, etc.)

###  File-System Under Linux
<p align="center">
<img src="https://github.com/Mo7ammedFarahat/MASRI-Jan25/blob/main/images/FileSys.png?raw=true" width="650"/>
</p>

---

<p align="center">
<img src="https://github.com/Mo7ammedFarahat/MASRI-Jan25/blob/main/images/FileSys2.png?raw=true" width="650"/>
</p>

### Home Directory and Working Directory
- When you first log in on a UNIX system, the working directory is your home directory.
- While working, you will be associated with one directory called the working directory or the current directory.
- An abbreviation of the working directory is displayed as part of the prompt on your terminal.
- The command `pwd` gives the absolute path of the working directory.

### What is a Path or a Pathname?

- A path locates a given file in the system hierarchy.
- An **absolute path** in the file system hierarchy for a given file or folder describes the parents all the way up to the root.
- A **relative path** describes the path to the file starting from the **current working directory**.

### ~ (Your Home Directory)

- `~` refers to the home directory in a given file system.
- The tilde `~` character can be used to specify paths starting at your home directory.

<p align="center">
<img src="https://github.com/Mo7ammedFarahat/MASRI-Jan25/blob/main/images/FileSys3.png?raw=true" width="650"/>
</p>

---

<p align="center">
<img src="https://github.com/Mo7ammedFarahat/MASRI-Jan25/blob/main/images/FileSys4.png?raw=true" width="650"/>
</p>

## Refer to the Parent and Current Directories

Every directory has two special subdirectories:

 `.` (dot): Refers to the current directory.  
 `..` (dot-dot): Refers to the parent directory.


<p align="center">
<img src="https://github.com/Mo7ammedFarahat/MASRI-Jan25/blob/main/images/FileSys5.png?raw=true" width="350"/>
</p>

---
<p align="center">
<img src="https://github.com/Mo7ammedFarahat/MASRI-Jan25/blob/main/images/FileSys6.png?raw=true" width="650"/>
</p>

---
## Creating directories and navigating

### Commands for Manipulating Directories

| Command | Description                                           |
|---------|-------------------------------------------------------|
| `mkdir` | Make directory: creates a new directory              |
| `rmdir` | Removes a directory                                  |
| `pwd`   | Displays the absolute path of the current working directory |
| `cd`    | Change directory: allows moving from one directory to another |
| `ls`    | Lists a directory's content                          |

### pwd Command

- `pwd`: Print Working Directory.
- Displays the absolute path of your current location in the file system.
- Try `pwd` on your terminal.
- Example output: `/home/YourUsername`

### ls Command

- `ls` lists the content of the current directory by default.
- Command structure: `ls [OPTION] [dirname]`.

##### Useful Options:
- `-l`: Shows sizes, modified date and time, file or folder name, owner of file, and permissions.
- `-a`: Lists all files, including hidden files starting with `.`.
- `-lh`: Shows sizes in an easily readable format.
- `-R`: Recursively lists sub-directories.
- `-lS`: Sorts by file sizes.

### Creating a Directory

- `mkdir`: Makes a directory.
- Command structure: `mkdir dirname [path]`.
- `mkdir dirname`: Creates a directory with the specified `dirname` in your current working directory.
- To create a directory elsewhere, specify the path: `mkdir dirname path`.

### Commands Basic Structure

The general structure of commands in Unix/Linux is:

#### Examples:
- `ls –lh /home/Watson/NU`: Lists the contents of `/home/Watson/NU` in a human-readable format.
- `pwd`: Prints the absolute path of the current working directory.
- `mkdir Test1`: Creates a directory named `Test1` in the current working directory.

### What You Should Know About File Names in Linux

- There is no real distinction between the names of ordinary files and the names of directory files.
- No two files in the same directory can have the same name.
- Files in different directories can have the same name.
- Linux is case-sensitive: `Sanger`, `sanger`, and `SANGER` are different and represent three distinct files.
- In most cases, file extensions (e.g., `.txt`, `.exe`) are optional.

### Move in the File System

- `cd`: Changes the working directory.
- Command structure: `cd <path>`.
- Specify the path name of the directory you want to move to.
- The path can be specified as either:
  - **Absolute path**: The full path from the root directory.
  - **Relative path**: The path relative to your current working directory.


### Remove a Directory

- `rmdir`: Removes a directory.
- Command structure: `rmdir dirname [path]`.
- It removes the `dirname` directory.
- The directory should be in your current working directory.
- To remove it from elsewhere, specify the path: `rmdir dirname path_to_directory`.
- `rmdir` only works if the directory is empty.
- If the directory contains files or subdirectories, an error message will appear: `Directory not empty`.
- To remove a directory and its contents, use `rm -r` (recursive), which will recursively remove the directory and all its contents.

### How to Get Help for a Command from the Terminal

 `man commandname`: In Linux, this command is used to display the user manual for any command that can be run in the terminal.  

## Useful shortcuts
- `cd`: When used with no arguments, `cd` changes the working directory to your home directory.  
- `cd ~user_name`: Moves to the specified user's home directory.
- `Ctrl+A`: Move the cursor to the beginning of the command line.
- `Ctrl+C`: End a running program and return to the prompt.
- `Ctrl+D`: Logout from the current shell session; equivalent to `exit`.
- `Tab`: Autocomplete a file name.
- `Tab` + `Tab`: Displays command completion possibilities.
- `Ctrl+L`: Clear the terminal.

---

## Manipulating files, useful commands and tips

### Basics of Manipulating File Commands

#### `touch` Command
- `touch` is used to create, change, and modify timestamps of a file.
- The `touch` command creates an empty (zero-byte) new file.
- Command structure: `touch filename`.
- You can create more than one file at a time: `touch filename1 filename2 filename3`.

#### `touch` Command Options
- `-a`: Change the access time only.
- `-c`: If the file does not exist, do not create it.
- `-d`: Update the access and modification times.
- `-m`: Change the modification time only.
- `-r`: Use the access and modification times of another file.
- `-t`: Create a file using a specified time.

For more details, check this: [touch Command in Linux](https://phoenixnap.com/kb/touch-command-in-linux)

### Text Editors

- `nano`: A simple and easy-to-use text editor.
  - Installed by default in Ubuntu and many other Linux distributions.
  - It’s a WYSIWYG editor: "What you see is what you get." What you type directly goes into the text input.
- `vim`, `emacs`, `gedit`, `Geany`: Excellent programs, but they require some learning.

<p align="center">
<img src="https://github.com/Mo7ammedFarahat/MASRI-Jan25/blob/main/images/Nano.png?raw=true" width="650"/>
</p>

### Get Started with nano

1.  `nano file1`: Open the file `file1` in nano.  
2.  Type: `"my first test file with webminal"`.  
3.  Hit **Enter** to move to another line and type: `"the second line of test"`.  
4.  Once you finish typing, press **Ctrl+X** to exit.  
5.  When prompted to save the modified buffer, answering "No" will destroy the changes.  
     Press **Y** to save.  
6. `nano file2`: Open a new file `file2`.  
7.  Type: `"my second test file with webminal"` and add any other 4 lines of text. 
  
  
    

####
To search for a text string, hit **Ctrl+W**, and enter your search term.  
This search can then be cancelled mid-execution by hitting **Ctrl+C** without destroying your buffer.  
**Ctrl+X**: Finish typing and close an open file.

#### Remember:
- `nano pathname`: Opens the file if it already exists, allowing you to modify and save changes.
- If the file does not exist, it creates a new file in the specified path.

### Displaying Whole Content of a File or Parts of It

- `cat`: View the content of a short file.  
  Syntax: `cat <filename>`

- `more`: View the content of a long file and navigate through it.  
  Syntax: `more <filename>`

- `less`: View the content of a long file, by portions.  
  Syntax: `less <filename>`

- `head`: View the first lines of a long file.  
  Syntax: `head <filename>`

- `tail`: View the last lines of a long file.  
  Syntax: `tail <filename>`

### View File Content: `less` Command

- The `less` command displays a text file's content, one page at a time.
- Structure: `less <filename>`
- Move a page down: Either use the **Page Down** key or **Space**.
- To exit `less`, type **q**.
- To go to the end of the text file, type **g**.

### `head` and `tail` Commands

- The `head` command displays a text file's content, by default:  
  10 first lines at a time.  
  Syntax: `head <options> <filename>`

- The `tail` command displays a text file's content, by default:  
  10 last lines at a time.  
  Syntax: `tail <options> <filename>`

### Copy, Move, and Remove

- `cp`: Copy files and directories.  
  Structure: `cp <pathfrom> <pathto>`

- `mv`: Move or rename files and directories.  
  Structure: `mv <pathfrom> <pathto>`

- `rm`: Remove files and directories.  
  Structure: `rm <pathname>`

### Copying Command: `cp`

- Simplest form: `cp file1 file2`  
  ➔ Copies the contents of `file1` into `file2`. If `file2` does not exist, it is created. Otherwise, `file2` is silently overwritten with the contents of `file1`.

- `cp filename dirpath`  
  ➔ Makes a copy of the file (or directory) into the specified destination directory.

### Other Examples: `cp`

   - Add the interactive mode with the option `-i`:  
     `cp -i file1 file2`  
     ➔ Same as the previous command. However, if `file2` exists, the user is notified before overwriting `file2` with the contents of `file1`.

   - `cp -R pathdir1 pathdir2`  
     ➔ Copies the contents of the directory `dir1`. If directory `dir2` does not exist, it is created. Otherwise, it creates a directory named `dir1` within directory `dir2`.

### Moving Command: `mv`

The `mv` command moves or renames files and directories, depending on how it is used.

- **To rename a file:**  
  `mv filename1 filename2`  
  ➔ If `file2` exists, its contents are silently replaced with the contents of `file1`. To avoid overwriting, use the interactive mode:  
  `mv -i filename1 filename2`

- **To move a file (or a directory) to another directory:**  
  `mv file dirpath`

- **To move different files (or directories) to another directory:**  
  `mv file1 file2 file3 dirpath`

- **To move a directory to another directory:**  
  `mv dir1 dir2`  
  ➔ If `dir2` does not exist, then `dir1` is renamed to `dir2`. If `dir2` exists, the directory `dir1` is moved within directory `dir2`.

### The `rm` Command

The `rm` command deletes files and directories.

- **To remove a file:**  
  `rm filename`

- **To remove many files:**  
  `rm filename1 filename2`

- **Add the interactive mode to prompt the user before deleting with `-i`:**  
  `rm -i filename1 filename2`

- **Delete directories with all their contents:**  
  `rm -r dir1 dir2`

### :warning: Be Careful with `rm`! :warning:

- Linux does not have an undelete command.
- Once you delete something with `rm`, it's gone!
- You can inflict terrific damage on your system with `rm` if you are not careful, particularly with wildcards.
- **Try this trick before using `rm`:** Construct your command using `ls` first to see what files will be affected.

### Wildcards

- Since the shell uses filenames so much, it provides special characters to help rapidly specify groups of filenames.
- A group of special characters called **wildcards** allows selecting filenames based on a pattern of characters.

### Wildcard Meanings

| Wildcard             | Meaning                                                                                                         |
|----------------------|-----------------------------------------------------------------------------------------------------------------|
| `*`                  | Matches any characters                                                                                          |
| `?`                  | Matches any single character                                                                                     |
| `[!characters]`      | Matches any character that is not a member of the set characters                                                 |
| `[characters]`       | Matches any character that is a member of the set characters. The set of characters may also be expressed as a  |
|                      | POSIX character class such as one of the following:                                                             |
|                      | `[:alnum:]` - Alphanumeric characters                                                                             |
|                      | `[:alpha:]` - Alphabetic characters                                                                              |
|                      | `[:digit:]` - Numerals                                                                                           |
|                      | `[:upper:]` - Uppercase alphabetic characters                                                                   |
|                      | `[:lower:]` - Lowercase alphabetic characters                                                                   |
| `a*`                 | Any file name starting with "a"                                                                                 |
| `*`                  | All possible filenames                                                                                          |
| `A*.fasta`           | All filenames that begin with "A" and end with `.fasta`                                                         |
| `????.vcf`           | Any filenames that contain exactly 4 characters and end with `.vcf`                                              |
| `[abc]*`             | Any filename that begins with "a" or "b" or "c" followed by any other characters                                |
| `[[:upper:]]*`       | Any filename that begins with an uppercase letter. This is an example of a character class                     |

---
### Download Files from the Web


- `wget` stands for "web get". It is a command-line utility which downloads files over a network.
- It supports HTTP, HTTPS, and FTP protocols.

**Syntax**:  
`wget [–options] [URL]`

### Let's Try It:

To get the fasta file of *Plasmodium falciparum* from PlasmoDB, use the following command:

```bash
wget http://plasmodb.org/common/downloads/release-9.0/Pfalciparum/fasta/PlasmoDB-9.0_Pfalciparum_BarcodeIsolates.fasta
```

---
## Basic operations on files and data extraction
#### Some Statistics About Your File Content: `wc` Command

- `wc` prints newline, word, and byte counts for each file.

**Syntax**:  
`wc <options> <filename>`

### Some Useful Options:
- `-c`: Prints the byte counts.
- `-m`: Prints the character counts.
- `-l`: Prints the newline counts.

For more information about the `wc` command, use the following command in the terminal:  

```bash
man wc
```
### Basic Operations on Files

- **`sort`**: Reorder the content of a file "alphabetically".  
  **Syntax**:  
  `sort <filename>`

- **`uniq`**: Removes duplicated lines from a file.  
  **Syntax**:  
  `uniq <filename>`

- **`join`**: Compare the contents of two files and output the common entries.  
  **Syntax**:  
  `join <filename1> <filename2>`

- **`diff`**: Compare the contents of two files and output the differences.  
  **Syntax**:  
  `diff <filename1> <filename2>`
### Sorting Data

- **`sort`**: Outputs a sorted order of the file content based on a specified sort key (default: sorts the entire input).  
  **Syntax**:  
  `sort <options> <filename>`

- **Default field separator**: Blank space.

- Sorted files are often used as input for other commands, so `sort` is commonly used in combination with other commands.

- For more options, use `man sort` to see the manual for detailed usage.

### Sorting Data: Examples

- **Sort alphabetically (default option)**:  
  `sort <filename>`

- **Sort numerically**:  
  `sort -n <filename>`

- **Sort on a specific column (e.g., column number 4)**:  
  `sort -k 4 <filename>`

- **Sort based on a tab separator**:  
  `sort -t $'\t' <filename>`

### Extracting Data from Files

- **grep**: Used to search for the occurrence of a specific pattern (regular expression using wildcards) in a file.  
  **Syntax**: `grep <pattern> <filename>`

- **cut**: Used to extract specific fields from a file.  
  **Syntax**: `cut <options> <filename>`

### grep Command

- **grep** ("global regular expression profile") is used to search for the occurrence of a specific pattern (regular expression) in a file.
- `grep` outputs the entire line containing the matching pattern.

**Syntax**: `grep <pattern> <filename>`

#### Examples:
- Extract lines containing the pattern `xxx` from a file:  
  `grep xxx <filename>`

- Extract lines that do not contain the pattern `xxx` from a file:  
  `grep -v xxx <filename>`

### `grep` example

Let’s consider a file named `ghandi.txt`:

```bash
$ cat ghandi.txt
The difference between what we do
and what we are capable of doing
would suffice to solve
most of the world's problems
```
```bash
$ grep what ghandi.txt
The difference between what we do
and what we are capable of doing
```
```bash
$ grep -v what ghandi.txt
would suffice to solve
most of the world's problems
```

### `cut` command

- `cut` is used to extract specific fields from a file.
- Structure: `cut <options> <filename>`
- For `<options>`, see `man`.
- Important options:
  - `-d` (field delimiter)
  - `-f` (field specifier)

#### Example:
Extract fields 2 and 3 from a file having ‘space’ as a separator:

```bash
cut -d' ' -f2,3 <filename>
```
### `uniq` command

- `uniq` outputs a file with no duplicated lines.
- `uniq` requires a sorted file as an input.
- Syntax: `uniq <options> <sorted_filename>`
- For `<options>`, see `man`.
- Useful option:
  - `-c` to output each line with its number of repeats.

### `join` command

- `join` is used to compare 2 input files based on the entries in a common field (called the "join field") and outputs a merged file.
- `join` requires sorted files as input.
- Lines with identical "join field" values will appear only once in the output.
- Syntax: `join <options> <filename1> <filename2>`
- For `<options>`, see `man`.

### `diff` command

- `diff` is used to compare 2 input files and displays the different entries.
- It can be used to highlight differences between 2 versions of the same file.
- Default output: Common lines are not shown; only different lines are indicated, showing what has been added (a), deleted (d), or changed (c).
- Syntax: `diff <options> <filename1> <filename2>`
- For `<options>`, see `man`.
---

## Outputs redirection and combining different commands
### Commands Outputs

- By default, the standard output of any command will appear on the terminal screen.
- Redirection of the output result to a file is possible.
- This is particularly useful for big files.
- Syntax: `command options filename.in > filename.out`
### Output redirection
- If the file exists, the result will be redirected to it.  
- If the file does not exist, it will be automatically created and the result redirected to it  
```bash
$ cat ghandi.txt
The difference between what we do
and what we are capable of doing
would suffice to solve
most of the world's problems
```
```bash
$ cut -d’ ‘ -f2,3 ghandi.txt
difference between
what we
suffice to
of the
```
```bash
$ cut -d’ ‘ -f2,3 ghandi.txt > ghandi.txt.out
```
```bash
$ cat ghandi.txt.out
difference between
what we
suffice to
of the
```
### Commands Combination

- The standard output of any command will be one unique output.
- As seen previously, this output can be printed on the screen or redirected to a file.
- However, the output result of a command can also be redirected to another command.
- This is particularly useful when several operations are needed for a file, with no need to store the intermediate outputs.

### Commands Combination: Example

- Combining several commands is done thanks to the use of a “|” character.
- **Structure:**

```bash
command1 options1 filename1.in | command2 options2 > filename.out
```

- This can be done for as many commands as needed.

## Permissions, groups and control
### Files and directories permissions
On a Linux system, each file and directory is assigned access rights for the owner of the file, the members of a group of related users, and everybody else.

```bash
mohammedfarahat@slurm-login:~$ ls -l
drwxr-xr--  2 mohammedfarahat cbio-group           8 Mar 17  2024 MolDock
drwxr-xr-x  4 mohammedfarahat cbio-group           2 Dec 18  2023 Projects
drwxr-xr-x  3 mohammedfarahat cbio-group           2 Jan 26  2024 PurB
drwxr-xr-x  3 mohammedfarahat cbio-group           1 Jul 17  2023 R
drwxr-xr-x  6 mohammedfarahat cbio-group          15 Mar  5  2024 RNASeq_V0
drwxrwxr-x  8 mohammedfarahat cbio-group           7 Jan 23  2024 RefGraph
```
as you see above the first section is the file permissions, followed by the link count of the directory or file (including itself), then the owner, the group, size in bytes, modification time and the file name.  
#### Permissions are broken into 4 sections:
<p align="center">
<img src="https://github.com/Mo7ammedFarahat/MASRI-Jan25/blob/main/images/FilePermission.png?raw=true" width="650"/>
</p>

### Access Permissions on Files

- **r**: Indicates **read** permission. It allows the user to read or copy the file's contents.
- **w**: Indicates **write** permission. It allows the user to modify or change the file.
- **x**: Indicates **execute** permission. It allows the user to run the file as a program or script, where appropriate.

### Access Permissions on Directories

- **r**: Indicates permission to **list** files in the directory.
- **w**: Indicates permission to **delete** files from the directory or **move** files into it.
- **x**: Indicates permission to **access** files in the directory. This means you may read files in the directory, provided you have read permission on the individual files.

### `chmod` command

- The `chmod` command is used to **change the permissions** of a file or a directory.
- **Syntax**: `chmod options permissions filename`
- Only the **owner** of the file can use `chmod` to change the permissions.
- Permissions define access for:
  - **Owner** (user who owns the file)
  - **Group** (users who are members of the file's group)
  - **Others** (everyone else)

#### Two Ways to Specify Permissions:
- **Symbols**: Alphanumeric characters
- **Octals**: Digits (0 to 7)

### `chmod` Options

| Symbol | Meaning                           |
|--------|-----------------------------------|
| u      | user                              |
| g      | group                             |
| o      | other                             |
| a      | all                               |
| r      | read                              |
| w      | write (and delete)                |
| x      | execute (and access directory)    |
| +      | add permission                    |
| -      | take away permission              |

### Octal Permissions

- 4 stands for "read"
- 2 stands for "write"
- 1 stands for "execute"
- 0 stands for "no permission"
### chmod Examples

 `chmod u=rwx,g=rx,o=r filename`  
  This is the same as:  
 `chmod 754 filename`
<p align="center">
<img src="https://github.com/Mo7ammedFarahat/MASRI-Jan25/blob/main/images/FilePermission2.png?raw=true" width="350"/>
</p>

#### More Examples

- `777` : `(rwxrwxrwx)` No restrictions on permissions. Anybody may do anything
- `755` : `(rwxrxrx)` The file's owner may read, write, and execute the file. All others may read and execute the file (common for programs that are used by all users)
- `700` : `(rwx)` The file's owner has all the rights. Nobody else has any rights (private for the owner)
- `666` : `(rwrwrw)` All users may read and write the file.
- `644` : `(rwrr)` The owner may read and write a file, while all others may only read the file (everybody may read, but only the owner may change)
- `600` : `(rw)` The owner may read and write a file. All others have no rights

## Environment variable
### Variables

- Variables are areas of memory that can be used to store information and are referred to by a name.
- **How to create a variable**: a line that contains the name of the variable followed immediately by an equal sign (`=`).
- 2 types of variables:
  - **Shell variables**
  - **Environment variables**
- Some variables are already set in your shell session.
- `printenv`: prints the values of all your environment variables.

### What is an Environment Variable

- An environment variable is a dynamic "object" on a computer that stores a value, which in turn can be referenced by one or more programs.
- Environment variables help programs know what directory to install files in, where to store temporary files, where to find user profile settings, and other things.
- Environment variables help to create and shape the environment where a program runs.
### Examples of Environment Variables

- **HOME**: The environmental variable that shows the current user's home directory.
- **PATH**: The environmental variable, which contains a colon-separated list of the directories that the system searches to find the executable program corresponding to a command issued by the user.
- **PWD**: Always stores the value of your current working directory.

## Shell Scripting
### echo Command

- **Syntax**: `echo options arguments`
- Writes arguments to the standard output.
- `echo`: Just prints its command-line parameters to standard output.
- If you redirect the result, your arguments will be written into the file you are redirecting to.
- Commonly used by shell scripts to display results or ask the user to enter parameters.
### Let’s Echo Some Stuff

- `echo Bioinformatics is great starting writing scripts`
- If you want to jump to another line, add `\n` and use the option `-e`:
  - `echo -e "Bioinformatics is great \n starting writing scripts"`
- Setting a variable: `X=firstvariable`
- `echo X`: prints `X`
- `echo $X`: prints `firstvariable` (the value of the variable)
- `echo '$X'`: ➔ `$X`
- `echo "$X"`: ➔ `firstvariable`

### Print the Result of a Command

- Asking the shell to substitute the results of a given command:
  - `` command `` or `$(command)`
- Example:
  - `echo \`pwd\`` or `echo $(pwd)`

### What is a Shell Script

- Short programs written in shell programming language useful to automate tasks under Linux OS.
- A shell script is a file containing a series of commands.
- Shell scripts could be helpful to perform the same actions on many different files.
- Shell script = scripting interpreter + command line interface to the system.
- `echo` is also commonly used in shell scripts to display a message or instructions, such as "Please enter Y or N" in an interactive session with users.

### Let’s Start Using the Power of Scripting

1. Create a new script file using `nano`:
```bash
   nano myfirstscript
```
2. Write the content of your script, for example:
```bash
clear
echo "Hey, I am starting writing shell scripts"
echo "Let the fun begin!!!"
```
3. Run your script using `./`:
```bash
./myfirstscript
```
4. Change the rights to make sure you have the right to execute:
```bash
chmod u+x myfirstscript
```
or
```bash
chmod 744 myfirstscript
```
### The Shebang

A shebang (`#!/bin/bash`) indicates the beginning of a script and specifies the program interpreter to use for running the script.

- The shebang is followed by the absolute path to the executable program.
- For a Perl script, the shebang could be `#!/usr/bin/perl`.
- You can use the `which` command to locate the executable file associated with a given command:
  ```bash
  which perl  # ➔ /usr/bin/perl
  which bash  # ➔ /bin/bash
  ```
### Use Variables in Your Scripts

- Variables make your script easier to maintain.
- They help reduce the amount of typing and make your script more flexible.

### If Statements in Shell Scripting

Syntax:

```bash
if [ conditional_expression ]
then
  commands
else
  commands
fi
```
Example:
```bash
#!/bin/bash
echo "Let’s try some conditional tests"
x=`find *.fasta | wc -l`
echo "The current working directory contains $x fasta files"
if [ $x -gt $y ]  # or if (($x > $y))
then
  echo "There are many existing fasta files in this directory"
else
  echo "There are very few fasta files in this directory. Here is the listing: `ls *.fasta`"
fi
```

### Loops in Shell Scripting (for)

#### Syntax:
```bash
for variable in values
do
  commands
done
```
```bash
#!/bin/bash
for x in file1 file2
do
  head -n 3 $x
  echo "Operation completed on file: $x"
done
```

### Loops in Shell Scripting (while)

#### Syntax:
```bash
while [ condition ]
do
  command1
  command2
done
```
```bash
#!/bin/bash
n=1
while [ $n -le 5 ] # n should have an initial value
do
  echo "Welcome $n times"
  n=$(( n + 1 )) # increment 
  $n
done
```

### Operators Supported by Shell

| Operator | Description                                                                                 | Example                                    |
|----------|---------------------------------------------------------------------------------------------|--------------------------------------------|
| `+`      | Addition - Adds values on either side of the operator.                                      | `expr $a + $b` will give 30.               |
| `-`      | Subtraction - Subtracts the right-hand operand from the left-hand operand.                  | `expr $a - $b` will give -10.              |
| `*`      | Multiplication - Multiplies values on either side of the operator.                          | `expr $a \* $b` will give 200.             |
| `/`      | Division - Divides the left-hand operand by the right-hand operand.                         | `expr $b / $a` will give 2.                |
| `=`      | Assignment - Assigns the right operand's value to the left operand.                         | `a=$b` would assign the value of `b` to `a`. |
| `==`     | Equality - Compares two numbers. If both are the same, returns true.                        | `[ $a == $b ]` would return false.         |
| `!=`     | Not Equality - Compares two numbers. If both are different, returns true.                   | `[ $a != $b ]` would return true.          |
| `-eq`    | Checks if the value of two operands are equal. If yes, the condition becomes true.          | `[ $a -eq $b ]` is not true.               |
| `-ne`    | Checks if the value of two operands are not equal. If not equal, the condition becomes true.| `[ $a -ne $b ]` is true.                   |
| `-gt`    | Checks if the value of the left operand is greater than the value of the right operand. If yes, the condition becomes true. | `[ $a -gt $b ]` is not true.  |
| `-lt`    | Checks if the value of the left operand is less than the value of the right operand. If yes, the condition becomes true.    | `[ $a -lt $b ]` is true.       |
| `-ge`    | Checks if the value of the left operand is greater than or equal to the value of the right operand. If yes, the condition becomes true. | `[ $a -ge $b ]` is not true. |
| `-le`    | Checks if the value of the left operand is less than or equal to the value of the right operand. If yes, the condition becomes true.    | `[ $a -le $b ]` is true.       |

## Controlling tasks
### Commands to Control Processes

- **ps**: List the processes running on the system.
- **kill**: Send a signal to one or more processes (usually to "kill" a process).
- **jobs**: An alternate way of listing your own processes.
- **bg**: Put a process in the background.

### Launching a Background Job

- Programs that take time or open a new Graphical User Interface may prevent the prompt from reappearing after the program is launched.  
  ➔ The shell waits for the program to finish before control returns to you.

#### Options to Interrupt or Background a Program

- **Ctrl+Z**: Interrupts a program.  
- **Run in the background**: Use the command name followed by `&` to put the program in the background so the prompt returns immediately.

## SSH into remote machine
#### What is SSH

**SSH (Secure Shell)** is a protocol used to securely log onto remote systems, such as remote Linux machines and Unix-like servers.

### SSH Command

The `ssh` command is the tool used in Linux to connect via the SSH protocol.

### Syntax

```bash
ssh remoteusername@remotehost
```
- Remote host could be an IP address or domain name
- You will be asked to provide your password
- To exit and go back to your into your local session, use `exit`.

### Multi-users in a Linux Machine

- While your computer only has one keyboard and monitor, it can still be used by more than one user.
- For example, if your computer is attached to a network or the Internet, remote users can log in via **SSH** (Secure Shell) and operate the computer.
- Remote users can execute applications and have the output displayed on a remote computer.

### Copy Files from or to a Remote Machine

- **scp**: Secure Copy
- **Syntax**: `scp pathfrom pathto`
- The difference: In **scp**, at least the source or the destination is on a remote machine.
- Example: Uploading all the `.txt` files from your current working directory to a remote host:
  ```bash
  scp ./*.txt username@myhost.com:/home/username/folder
  ```

