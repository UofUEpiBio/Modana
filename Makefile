build: NAMESPACE
	R CMD build . 

check: install
	R CMD check --no-manual *_*.tar.gz

install: build
	R CMD INSTALL *_*.tar.gz
	
  
NAMESPACE: R/*
	Rscript --vanilla -e 'roxygen2::roxygenize()'

docs: NAMESPACE
  
.PHONY: build check install docs

covr: 
	rm -rf covr || mkdir covr && \
		Rscript -e 'pth <- normalizePath("~/Modana/covr");Sys.setenv(SLURMR_TMP_PATH=pth);saveRDS(covr::package_coverage(quiet = FALSE, clean = FALSE, install_path = pth), "covr/dat.rds")' && \
		echo "Now you can go ahead and upload the coverage"

