## This is MeanRarity

current: target
-include target.mk

# -include makestuff/perl.def

vim_session:
	bash -cl "vmt"

######################################################################

Sources += .Rbuildignore

Ignore += Meta/ doc/

package:
	R CMD INSTALL .

######################################################################

Sources += $(wildcard dev/*.R)

dev/updown.Rout: dev/updown.R
	$(pipeR)

dev/udtest.Rout: dev/udtest.R dev/updown.rda
	$(pipeR)

######################################################################

## pages/ currently

######################################################################

Sources += vignettes/*.Rmd

%.html: %.Rmd
	$(rmdh_r)

## vignettes/Using_MeanRarity.vig.html: vignettes/Using_MeanRarity.Rmd
vignettes/%.vig.html: vignettes/%.Rmd
	$(rmdh_r)

######################################################################

### Makestuff

Sources += Makefile

Ignore += makestuff
msrepo = https://github.com/dushoff

## Want to chain and make makestuff if it doesn't exist
## Compress this ¶ to choose default makestuff route
Makefile: makestuff/Makefile
makestuff/Makefile:
	git clone $(msrepo)/makestuff
	ls makestuff/Makefile

-include makestuff/os.mk

-include makestuff/pipeR.mk
-include makestuff/rpkg.mk
-include makestuff/rmd.mk

-include makestuff/git.mk
-include makestuff/visual.mk
