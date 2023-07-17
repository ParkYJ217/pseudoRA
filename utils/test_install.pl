#!/usr/bin/env perl

use strict;
use warnings;

my $scriptdir;
my $ecode;


# Where are we running from?
use File::Basename;
use Cwd 'abs_path';

our $warnings;

print("\n");
print("Checking the OS\n");
my $os = "$^O";
if($os eq "linux") { print("\t$os OK\n"); } else { die("The \"$os\" OS was detected, but PseudoRA only runs on Linux!!\n\n"); }


print("\n");
print("Checking that samtools is installed\n");
check_command("samtools --help", "ERROR: The samtools program can not be found in this environment!!");


print("\n");
print("Checking that bwa is installed\n");
check_command("man bwa", "ERROR: The bwa interpreter can not be found in this environment!!");


print("\n");
print("Checking that java is installed\n");
check_command("java -h", "ERROR: The java VM can not be found in this environment!!");


print("\n");
print("Checking that all the required python libraries are available in this environment\n");
$ecode = check_command("python3 -h", "The python3 interpreter can not be found in this environment!!");
if(!$ecode) {
	check_python_library("numpy");
	check_python_library("pandas");
	check_python_library("pysam");
	check_python_library("Bio");
}


if($warnings) {
	print("\n");
	print("------------------------------------------------------------------------------\n");
	print("\n");
	print("WARNING: Some PseudoRA dependencies could not be found in your environment!\n");
	print($warnings);
	die("\n");
} else {
	print("\n");
	print("All checks successful\n");
}


print("\n");


sub check_command {
	my $command = $_[0];
	my $msg = $_[1];
	my $ecode = system("$command > /dev/null 2>&1");
        if(!$ecode) { print("\t$command OK\n"    ); }
	else {
		warn("\t$command NOT OK\n"); 
		if($msg) { print("\t\t$msg\n"); }
		$warnings .= "\t- $msg\n";
	}
	return $ecode;
}


sub check_python_library {
	my $library = $_[0];
	my $msg = $_[1];
	if(!$msg) { $msg = "Missing python library \"$library\""; }
	my $command = "python3 -c 'import $library'";
        return check_command($command, $msg);
}


