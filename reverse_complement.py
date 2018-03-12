from sys import argv
import re

script, text = argv

#print DNA sequence

print text

rep = {"A": "T", "T": "A", "G": "C", "C": "G"} # define desired replacements here

# use these three lines to do the replacement
rep = dict((re.escape(k), v) for k, v in rep.iteritems())
pattern = re.compile("|".join(rep.keys()))
rev_text = pattern.sub(lambda m: rep[re.escape(m.group(0))], text)


print rev_text[::-1]



