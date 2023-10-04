project: fenvelopes
summary: Calculation of phase envelops using Equations of State
project_github: https://github.com/fedebenelli/fenvelopes
author: Federico Benelli
author_description: PhD student with focus on reservoir PVT simulation.
author_email: federico.benelli@mi.unc.edu.ar
github: https://github.com/fedebenelli
src_dir: ../src 
         ../app
exclude_dir: ../test ../doc
output_dir: ../doc/ford_site
preprocessor: gfortran -E
display: public
         protected
         private
source: false
proc_internals: true
sort: permission-alpha
docmark_alt: -!
docmark: !
predocmark_alt: *
print_creation_date: true
creation_date: %Y-%m-%d %H:%M %z
md_extensions: markdown.extensions.toc
               markdown.extensions.smarty
graph: true
license: MIT

{!../README.md!}
