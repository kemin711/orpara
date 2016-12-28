#!/usr/bin/perl -w

# this program will eliminate the problem cause by
# ^M symbol appearing in your vi terminal for text files created on
# window's platform.
# useage unixnewline [input text file] > [output file]
while (<>) {
   s/\r/\n/g;
   print;
}

