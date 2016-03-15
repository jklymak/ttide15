nbs = $(wildcard *.ipynb)
htmls = $(patsubst %.ipynb,%.html,$(nbs))
pys = $(wildcard *.py)
ms = $(wildcard *.m)

all: index html jmkpython.tgz

index: index.html

index.html: README.md
	markdown README.md > index.html

html: $(htmls)

%.html: %.ipynb 
	jupyter-nbconvert $< --to html

jmkpython.tgz: $(pys) $(ms)
	tar cvzf jmkpython.tgz *.py *.m

install:
	rsync -av ~/ttide/ttide15/ ~/Sites/ttide15 --exclude .git
