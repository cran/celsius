all :: man check install
test :: check
man ::
	/usr/bin/perl ./man.pl
check ::
	/home/mcarlson/arch/i686/R-patched/bin/R CMD check -l /home/mcarlson/R/i686-pc-linux-gnu-library/2.5 `pwd`
install ::
	/home/mcarlson/arch/i686/R-patched/bin/R CMD INSTALL -l /home/mcarlson/R/i686-pc-linux-gnu-library/2.5 `pwd`
dist ::
	rm -f ./*.gz
	rm -f ./*.zip
	rm -rf ./*.Rcheck
	/home/mcarlson/arch/i686/R-patched/bin/R CMD build `pwd`
	/home/mcarlson/arch/i686/R-patched/bin/R CMD INSTALL -l . ./*.gz
	find `cat DESCRIPTION | grep -E '^Pack' | awk '{print $$2}'` | zip -@ `cat DESCRIPTION | perl -e '@F=<>;chomp(@F);%F=map{/^(\S+):\s+(\S+)$$/;$$1=>$$2}@F;print "$$F{Package}_$$F{Version}.zip"'`
	rm -rf `cat DESCRIPTION | grep -E '^Pack' | awk '{print $$2}'`
