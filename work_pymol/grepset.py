from pymol import cmd
import pymol.setting
import sys

def grepset(regexp='',output=sys.stdout):
   '''
DESCRIPTION
   "grepset" greps through the list of settings using a python
   regular expression as defined in the 're' module.
   It returns a list of settings/values matching the regexp.
   No regexp returns every setting.

USAGE
   grepset [regexp]

EXAMPLE
   grepset line
   grepset ray

SEE ALSO
   Python re module
   '''

   from re import compile

   if output != sys.stdout:
     out_handle=open(output,'w')
   else:
     out_handle=sys.stdout
   count=0
   regexp=compile(regexp)
   for a in pymol.setting.get_index_list():
      setting=pymol.setting._get_name(a)
      if regexp.search(setting):
         count = count + 1
         out_handle.write('%-30s %s\n' % (setting, cmd.get_setting_text(a,'',-1)))

   out_handle.write('%d settings matched\n' % count)
cmd.extend('grepset',grepset)
