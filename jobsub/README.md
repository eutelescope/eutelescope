jobsub is a tool for the convenient run-specific modification of
Marlin steering files and their execution through the Marlin processor.

Usage
===============================================================================
```
usage: jobsub [-h] [--option NAME=VALUE] [-c FILE] [-csv FILE] [-g]
              [-condor FILE] [-lx FILE] [--concatenate] [--log-file FILE]
              [-l LEVEL] [-s] [--dry-run] [--plain]
              jobtask [runs [runs ...]]

A tool for the convenient run-specific modification of Marlin steering files
and their execution through the Marlin processor

positional arguments:
  jobtask               Which task to submit (e.g. converter, hitmaker, align);
                        task names are arbitrary and can be set up by the
                        user; they determine e.g. the config section and
                        default steering file names.
  runs                  The runs to be analyzed; can be a list of single runs
                        and/or a range, e.g. 1056-1060.

optional arguments:
  -h, --help            show this help message and exit
  --option NAME=VALUE, -o NAME=VALUE
                        Specify further options such as 'beamenergy=5.3'. This
                        switch be specified several times for multiple options
                        or can parse a comma-separated list of options. This
                        switch overrides any config file options.
  -c FILE, --conf-file FILE, --config FILE
                        Load specified config file with global and task
                        specific variables
  -csv FILE, --csv-file FILE
                        Load additional run-specific variables from table
                        (text file in csv format)
  -g, --graphic
  -condor FILE, --condor_file FILE
                        Specify parameter file for HTCondor submission. Run
                        batch submission via condor_submit instead of calling
                        Marlin directly
  -lx FILE, --lxplus-file FILE, --lxplus FILE
                        Specify bsub parameter file for LXPLUS submission. Run
                        LXPLUS submission via bsub instead of calling Marlin
                        directly
  --concatenate         Modifies run range treatment: concatenate all runs
                        into first run (e.g. to combine runs for alignment) by
                        combining every options that includes the string
                        '@RunRange@' multiple times, once for each run of the
                        range specified.
  --log-file FILE       Save submission log to specified file
  -l LEVEL, --log LEVEL
                        Sets the verbosity of log messages during job
                        submission where LEVEL is either debug, info, warning
                        or error
  -s, --silent          Suppress non-error (stdout) Marlin output to console
  --dry-run             Write steering files but skip actual Marlin execution
  --plain               Output written to stdout/stderr and log file in
                        prefix-less format i.e. without time stamping
```

Preparation of Steering File Templates
===============================================================================
  Steering file templates are valid Marlin steering files (in ```xml```
  format) where single values are replaced by variables in the form
  ```@SomeVariable@```.

  When jobsub is run, these placeholders are filled with a
  user-defined value that can be specified through any of these
  sources (in order of precedence): command-line arguments, a config
  file, or a table with a row for each run number processed.
  
  There is only one predefined placeholder, @RunNumber@, which will be
  substituted with the current run number (padded with leading zeros
  to six digits, e.g. 001234).

Configuration
===============================================================================
  There are only very few predefined options: 
  * TemplateFile, TemplatePath: used to find the correct steering file template for the current task
  * DatabasePath, LcioPath, HistogramPath, LogPath, SteeringPath: paths for storing of results
  
  You can modify these options is the same way as placeholders in the template file, as described below.
  
  * Command Line
  
    Variable substitutions can be specified using the ```--option``` or ```-o``` command line switches, e.g.
    ```
    jobsub.py --option beamenergy=5.3 align 1234
    ```
    This switch be specified several times for multiple options or can
    parse a comma-separated list of options. This switch overrides any
    config file options.
   
   * Config File
   
     Config files are text file consisting of sections (indicated by '[]'):
      * a global section called ```[DEFAULT]```
      * task-specific sections
      * "name: value" or "name=value" entries, where 'name' are
      arbitrary steering file variables (case-insensitive).

     Some noteworthy features include:
      * comment prefix characters are # and ;
      * interpolation of format strings is supported, for example:
      ```
      [My Section]
      foodir: %(dir)s/whatever
      dir=frob
      long: this value continues
         in the next line
      ```
      would resolve the ```%(dir)s``` to the value of dir (frob in this case).
      * some default interpolations are ```%(home)s``` and ```%(eutelescopePath)s```
        which are set up with the environment variables ```$HOME``` and
        ```$EUTELESCOPE```, respectively.
      * The string ```@RunNumber@``` will be replaced in the template *after*
        all other variable strings were filled-in; therefore, you can use
        the ```@RunNumber@``` placeholder inside options, e.g. the file name.
        It will be replaced by the run number padded with leading zeros to 6 digits.
      * for more details, see the documentation to the Python module used
        for parsing: http://docs.python.org/2/library/configparser.html
      * for an example configuration file, please have a look on the provided examples
 
  * Table (comma-separated text file)
  
    * format: e.g. export from Open/LibreOffice with default settings (UTF-8,comma-separated, text-field delimiter: ") or emacs org-mode table (see http://orgmode.org/manual/Tables.html)
    * commented lines (starting with #) are ignored
    * first row (after comments) has to provide column headers which identify the variables in the steering template to replace (case-insensitive)
    * requires one column labeled "RunNumber"
    * only considers placeholders left in the steering template after processing command-line arguments and config file options
    
Concatenation
===============================================================================
If you have an option e.g. the LCIO input files that you want to
   fill with several runs in one steering file, you can use a command
   line switch to activate concatenation. This replaces any steering
   file placeholder whose corresponding option contains the string
   ```@RunRange@``` multiple times, once for every run specified.
   ```
   jobsub.py --concatenate --option LCIOInputFiles=/my/path/to/data/@RunRange@.lcio align 1234 1235-1237
   ```
   This will create *one* steering file (for run 1234) in which the placeholder
   ```@LCIOInputFiles@``` is replaced four times by its value with
   ```@RunRange@``` replaced by values from 1234 to 1237, e.g.
  ```
    <parameter name="FileName" type="string" value= @LCIOInputFiles@/>
   ```
   becomes:
   ```
    <parameter name="FileName" type="string" value= /my/path/to/data/1234.lcio 
            /my/path/to/data/1235.lcio /my/path/to/data/1236.lcio /my/path/to/data/1237.lcio/>
   ```
   This can be useful if you want to combine several runs e.g. for alignment.


Workflow
===============================================================================
  An analysis is controlled by a config file (config.cfg), a csv-table 
  (runlist.csv) and steering file templates (*.xml).

  In principle, neither the table nor the config are required as long
  as the template files do not contain any variables except for the
  run number (```@RunNumber@```).
  
  Execution step:
  ```
  jobsub -c config.cfg -csv runlist.csv JOBTASK RUNNR
  ```
  Here, jobsub will generate a steering file using the template file
  specified in the config file (default would be 'JOBTASK-tmp.xml'),
  thereby replacing any variables given in the config and table files.
  The final steering file will be processed by executing Marlin.
  
  The output paths of the analysis can be set inside the config file.
