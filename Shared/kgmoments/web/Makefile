all: codes web

.phony: codes web data

kgmoments.tar.gz : ../.git/index	
	(cd ..; git archive --prefix kgmoments/ master) | gzip > kgmoments.tar.gz

codes : kgmoments.tar.gz

index.md : ../README.md
	cp ../README.md index.md

web: index.md

deploy: all
	rsync -rtvp ./ ~/Dropbox/websites/purdue/codes/kgmoments/ --exclude=Makefile --exclude=*~ --exclude=.DS_Store

