doxygen = find_program('Doxygen')

if(doxygen.found())
  src_doxygen = [meson.source_root()+'/README.md', \
                 meson.source_root()+'/stable_distribution/adaptive_integration.h', \
                 meson.source_root()+'/stable_distribution/cubicspline.h', \
                 meson.source_root()+'/stable_distribution/gauss_kronrod.h', \
                 meson.source_root()+'/stable_distribution/myFloat.h', \
                 meson.source_root()+'/stable_distribution/stable_distribution.h', \
                 meson.source_root()+'/stable_distribution/stable_distribution_Vec.h', \
                 meson.source_root()+'/stable_distribution/stable_distribution_fit.h']

	doc_config = configuration_data()
	doc_config.set('PACKAGE_NAME', meson.project_name())
	doc_config.set('PACKAGE_VERSION', meson.project_version())
	doc_config.set('PROJECT_SOURCE_DIR', meson.source_root())
	doc_config.set('INPUT', ' '.join(src_doxygen))

	doxyfile = configure_file(input : 'Doxyfile.in',
				  output : 'Doxyfile',
				  configuration : doc_config,
				  install : false)
	custom_target('doxygen',
		      input : doxyfile,
		      output : ['html', 'latex']  ,
		      command : [ doxygen, doxyfile ],
		      depend_files: src_doxygen,
		      install : false,
		      build_by_default : true)

endif
