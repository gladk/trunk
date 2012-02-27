#!/usr/bin/scons
# coding: UTF-8
# vim:syntax=python:

#
# This is the master build file for scons (http://www.scons.org). It is experimental, though it build very well for me. Prequisities for running are having scons installed (debian & family: package scons)
#
# Type "scons -h" for yade-specific options and "scons -H" for scons' options. Note that yade options will be remembered (saved in scons.config) so that you need to specify them only for the first time. Like this, for example:
#
#	scons -j2 brief=1 debug=0 optimize=1 profile=1 exclude=extra,lattice,snow
#
# Next time, you can simply run "scons" or "scons -j4" (for 4-parallel builds) to rebuild targets that need it. IF YOU NEED TO READ CODE IN THIS FILE, SOMETHING IS BROKEN AND YOU SHOULD REALLY TELL ME.
#
# Scons will do preparatory steps (generating SConscript files from .pro files, creating local include directory and symlinking all headers from there, ...), compile files and install them.
#
# To clean the build, run `scons -c'. Please note that it will also _uninstall_ yade from $PREFIX!
#
# TODO:
#  1. [REMOVED] [DONE] retrieve target list and append install targets dynamically;
#  2. [DONE] configuration and option handling.
#  3. Have build.log with commands and all output...
#
# And as usually, clean up the code, get rid of workarounds and hacks etc.
import os,os.path,string,re,sys,shutil
import SCons
# SCons version numbers are needed a few times
# rewritten for python2.4
ver=[str(x) for x in SCons.__version__.split('.')]
for i in range(0,len(ver)):
	def any2num(x):
		if x in ('0123456789'): return x
		else: return '0'
	ver[i]="".join([any2num(x) for x in ver[i]])
	if len(ver[i])>2: ver[i]=ver[i][0:2]+"."+ver[i][2:]
sconsVersion=10000*float(ver[0])
if len(ver)>1: sconsVersion+=100*float(ver[1])
if len(ver)>2: sconsVersion+=float(ver[2])
# print sconsVersion
#

##########################################################################################
########## PROXY TO NEWER SCONS (DOWNLOADED IF NEEDED) ###################################
##########################################################################################
#print sconsVersion
if sconsVersion<9806.0 and not os.environ.has_key('NO_SCONS_GET_RECENT'):
	# sconsVersion<10200.0 
	# tgzParams=("http://dfn.dl.sourceforge.net/sourceforge/scons/scons-local-1.2.0.d20090223.tar.gz","/scons-local-1.2.0.d20090223")
	tgzParams=("http://heanet.dl.sourceforge.net/sourceforge/scons/scons-local-1.0.0.tar.gz","/scons-local-1.0.0")
	newPrefix="./scons-local";
	newUrl,newDir=tgzParams[0],newPrefix+"/"+tgzParams[1]
	if not os.path.exists(newDir):
		print "Scons version too old, downloading new version. All subsequent calls will be proxied to the new version transparently."
		import urllib,tarfile
		(filename,headers)=urllib.urlretrieve(newUrl)
		print filename,"\n",headers
		tar=tarfile.open(filename, "r:gz")
		for tarinfo in tar: tar.extract(tarinfo,newPrefix)
		print "Done extracting scons to",newDir
		assert(os.path.exists(newDir))
	if os.path.exists(newDir):
		Exit(os.execv(newPrefix+"/scons.py",[newPrefix+'/scons.py']+sys.argv[1:]))

##########################################################################################
############# OPTIONS ####################################################################
##########################################################################################

### hardy compatibility (may be removed at some later point probably)
### used only when building debian package, since at that point auto download of newer scons is disabled
if 'Variables' not in dir():
	Variables=Options
	BoolVariable,ListVariable,EnumVariable=BoolOption,ListOption,EnumOption

env=Environment(tools=['default','scanreplace'],toolpath=['scripts'])
profileFile='scons.current-profile'
env['sourceRoot']=os.getcwd()

profOpts=Variables(profileFile)
profOpts.Add(('profile','Config profile to use (predefined: default or ""); append ! to use it but not save for next build (in scons.current-profile)','default'))
profOpts.Update(env)
# multiple profiles - run them all at the same time
# take care not to save current profile for those parallel builds
if env['profile']=='': env['profile']='default'
# save the profile only if the last char is not !
saveProfile=True
if env['profile'][-1]=='!': env['profile'],noSaveProfile=env['profile'][:-1],False
if saveProfile: profOpts.Save(profileFile,env)

if env['profile']=='': env['profile']='default'
optsFile='scons.profile-'+env['profile']
profile=env['profile']
print '@@@ Using profile',profile,'('+optsFile+') @@@'

opts=Variables(optsFile)
## compatibility hack again
if 'AddVariables' not in dir(opts): opts.AddVariables=opts.AddOptions

def colonSplit(x): return x.split(':')

#
# The convention now is, that
#  1. CAPITALIZED options are
#   (a) notorious shell variables (they correspons most the time to SCons environment options of the same name - like CPPPATH)
#   (b) c preprocessor macros available to the program source (like PREFIX and SUFFIX)
#  2. lowercase options influence the building process, compiler options and the like.
#


opts.AddVariables(
	### OLD: use PathOption with PathOption.PathIsDirCreate, but that doesn't exist in 0.96.1!
	('PREFIX','Install path prefix','/usr/local'),
	('runtimePREFIX','Runtime path prefix; DO NOT USE, inteded for packaging only.',None),
	('variant','Build variant, will be suffixed to all files, along with version.','' if profile=='default' else '-'+profile,None,lambda x:x),
	BoolVariable('debug', 'Enable debugging information',0),
	BoolVariable('gprof','Enable profiling information for gprof',0),
	('optimize','Turn on optimizations (-1, 0 or 1); negative value sets optimization based on debugging: not optimize with debugging and vice versa.',-1,None,int),
	ListVariable('exclude','Yade components that will not be built','none',names=['gui','extra','common','dem','lattice','snow']),
	EnumVariable('PGO','Whether to "gen"erate or "use" Profile-Guided Optimization','',['','gen','use'],{'no':'','0':'','false':''},1),
	ListVariable('features','Optional features that are turned on','log4cxx,opengl,gts,openmp,vtk',names=['opengl','log4cxx','cgal','openmp','gts','vtk','python','gl2ps','devirt-functors','qt4','never_use_this_one']),
	('jobs','Number of jobs to run at the same time (same as -j, but saved)',2,None,int),
	#('extraModules', 'Extra directories with their own SConscript files (must be in-tree) (whitespace separated)',None,None,Split),
	('buildPrefix','Where to create build-[version][variant] directory for intermediary files','..'),
	('hotPlugins','Files (without the .cpp extension) that will be compiled separately even in the monolithic build (use for those that you modify frequently); comma-separated.',''),
	('chunkSize','Maximum files to compile in one translation unit when building plugins. (unlimited if <= 0, per-file linkage is used if 1)',7,None,int),
	('version','Yade version (if not specified, guess will be attempted)',None),
	('realVersion','Revision (usually bzr revision); guessed automatically unless specified',None),
	('CPPPATH', 'Additional paths for the C preprocessor (colon-separated)','/usr/include/vtk-5.0:/usr/include/vtk-5.2:/usr/include/vtk-5.4:/usr/include/vtk-5.6:/usr/include/eigen2:/usr/include/vtk'), # hardy has vtk-5.0
	('LIBPATH','Additional paths for the linker (colon-separated)',None),
	('libstdcxx','Specify libstdc++ location by hand (opened dynamically at startup), usually not needed',None),
	('QT4CXX','Specify a different compiler for files including qt4; this is necessary for older qt version (<=4.7) which don\'t compile with clang',None),
	('QT4DIR','Directory where Qt4 is installed','/usr/share/qt4'),
	('PATH','Path (not imported automatically from the shell) (colon-separated)',None),
	('CXX','The c++ compiler','g++'),
	('CXXFLAGS','Additional compiler flags for compilation (like -march=core2).',None,None,Split),
	('march','Architecture to use with -march=... when optimizing','native',None,None),
	BoolVariable('mono','[experimental] Build only one shared library and make all other files (python objects, for instance) only be symlinks.',0),
	('execCheck','Name of the main script that should be installed; if the current one differs, an erro is raised (do not use directly, only intended for --rebuild',None),
	('defThreads','Default number of OpenMP threads to run with, unless overridden with -j. Keep unset (negative) to default to using all cores.',-1),
	#('SHLINK','Linker for shared objects','g++'),
	('SHCCFLAGS','Additional compiler flags for linking (for plugins).',None,None,Split),
	BoolVariable('QUAD_PRECISION','typedef Real as long double (=quad)',0),
	BoolVariable('brief',"Don't show commands being run, only what files are being compiled/linked/installed",True),
)
opts.Update(env)

if str(env['features'])=='all':
	print 'ERROR: using "features=all" is illegal, since it breaks feature detection at runtime (SCons limitation). Write out all features separated by commas instead. Sorry.'
	Exit(1)

if saveProfile: opts.Save(optsFile,env)
# fix expansion in python substitution by assigning the right value if not specified
if not env.has_key('runtimePREFIX') or not env['runtimePREFIX']: env['runtimePREFIX']=env['PREFIX']
# set optimization based on debug, if required
if env['optimize']<0: env['optimize']=not env['debug']

# handle colon-separated lists:
for k in ('CPPPATH','LIBPATH','QTDIR','PATH'):
	if env.has_key(k):
		env[k]=colonSplit(env[k])

# do not propagate PATH from outside, to ensure identical builds on different machines
#env.Append(ENV={'PATH':['/usr/local/bin','/bin','/usr/bin']})
# ccache needs $HOME to be set, also CCACHE_PREFIX if used; colorgcc needs $TERM; distcc wants DISTCC_HOSTS
# fakeroot needs FAKEROOTKEY and LD_PRELOAD
propagatedEnvVars=['HOME','TERM','DISTCC_HOSTS','LD_PRELOAD','FAKEROOTKEY','LD_LIBRARY_PATH','CCACHE_PREFIX']
for v in propagatedEnvVars:
	if os.environ.has_key(v): env.Append(ENV={v:os.environ[v]})
if env.has_key('PATH'): env.Append(ENV={'PATH':env['PATH']})

# get number of jobs from DEB_BUILD_OPTIONS if defined
# see http://www.de.debian.org/doc/debian-policy/ch-source.html#s-debianrules-options
if os.environ.has_key('DEB_BUILD_OPTIONS'):
	for opt in os.environ['DEB_BUILD_OPTIONS'].split():
		if opt.startswith('parallel='):
			j=opt.split('=')[1].split(',')[0] # there was a case of 5, in hardy...
			print "Setting number of jobs (using DEB_BUILD_OPTIONS) to `%s'"%j
			env['jobs']=int(j)



if sconsVersion>9700: opts.FormatOptionHelpText=lambda env,opt,help,default,actual,alias: "%10s: %5s [%s] (%s)\n"%(opt,actual,default,help)
else: opts.FormatOptionHelpText=lambda env,opt,help,default,actual: "%10s: %5s [%s] (%s)\n"%(opt,actual,default,help)
Help(opts.GenerateHelpText(env))

###########################################
################# BUILD DIRECTORY #########
###########################################
import yadeSCons
## ALL generated stuff should go here - therefore we must determine it very early!!
if not env.has_key('realVersion') or not env['realVersion']: env['realVersion']=yadeSCons.getRealVersion() or 'unknown' # unknown if nothing returned
if not env.has_key('version'): env['version']=env['realVersion']

env['SUFFIX']='-'+env['version']+env['variant']
env['SUFFIX_DBG']=env['SUFFIX']+('' if not env['debug'] else '/dbg')
env['LIBDIR']='$PREFIX/lib/yade$SUFFIX_DBG'
print "Yade version is `%s' (%s), installed files will be suffixed with `%s'."%(env['version'],env['realVersion'],env['SUFFIX'])
buildDir=os.path.abspath(env.subst('$buildPrefix/build$SUFFIX_DBG'))
print "All intermediary files will be in `%s'."%env.subst(buildDir)
env['buildDir']=buildDir
# these MUST be first so that builddir's headers are read before any locally installed ones
buildInc='$buildDir/include'
env.Append(CPPPATH=[buildInc])

env.SConsignFile(buildDir+'/scons-signatures')

##########################################################################################
############# SHORTCUT TARGETS ###########################################################
##########################################################################################

if len(sys.argv)>1 and ('clean' in sys.argv) or ('tags' in sys.argv) or ('doc' in sys.argv):
	if 'clean' in sys.argv:
		if os.path.exists(buildDir):
			print "Cleaning: %s."%buildDir; shutil.rmtree(buildDir)
		else: print "Nothing to clean: %s."%buildDir
		sys.argv.remove('clean')
	if 'tags' in sys.argv:
		cmd="ctags -R --extra=+q --fields=+n --exclude='.*' --exclude=attic --exclude=doc --exclude=scons-local --exclude=include --exclude=lib/triangulation --exclude=debian --exclude='*.so' --exclude='*.s' --exclude='*.ii' --langmap=c++:+.inl,c++:+.tpp,c++:+.ipp"
		print cmd; os.system(cmd)
		sys.argv.remove('tags')
	if 'doc' in sys.argv:
		raise RuntimeError("'doc' argument not supported by scons now")
	#	cmd="cd doc; doxygen Doxyfile"
	#	print cmd; os.system(cmd)
	#	sys.argv.remove('doc')
	# still something on the command line? Ignore, but warn about that
	if len(sys.argv)>1: print "!!! WARNING: Shortcuts (clean,tags,doc) cannot be mixed with regular targets or options; ignoring:\n"+''.join(["!!!\t"+a+"\n" for a in sys.argv[1:]])
	# bail out
	Exit(0)


##########################################################################################
############# CONFIGURATION ##############################################################
##########################################################################################

# ensure non-None
env.Append(CPPPATH='',LIBPATH='',LIBS='',CXXFLAGS='',SHCCFLAGS='')

if 0:
	def CheckQt(context, qtdirs):
		"Attempts to localize qt3 installation in given qtdirs. Sets necessary variables if found and returns True; otherwise returns False."
		# make sure they exist and save them for restoring if a test fails
		origs={'LIBS':context.env['LIBS'],'LIBPATH':context.env['LIBPATH'],'CPPPATH':context.env['CPPPATH']}
		qtdirs=qtdirs[0].split()
		for qtdir in qtdirs:
			context.Message( 'Checking for qt-mt in '+qtdir+'... ' )
			context.env['QTDIR']=qtdir
			context.env.Append(LIBS='qt-mt',LIBPATH=qtdir+'/lib',CPPPATH=qtdir+'/include' )
			ret=context.TryLink('#include<qapplication.h>\nint main(int argc, char **argv){QApplication qapp(argc, argv);return 0;}\n','.cpp')
			context.Result(ret)
			if not ret:
				for k in origs.keys(): context.env[k]=origs[k]
			else:
				return ret
		return False


def CheckCXX(context):
	context.Message('Checking whether c++ compiler "%s" works...'%env['CXX'])
	ret=context.TryLink('#include<iostream>\nint main(int argc, char**argv){std::cerr<<std::endl;return 0;}\n','.cpp')
	context.Result(ret)
	return ret
def CheckLibStdCxx(context):
	context.Message('Finding libstdc++ library... ')
	if context.env.has_key('libstdcxx') and context.env['libstdcxx']:
		l=context.env['libstdcxx']
		context.Result(l+' (specified by the user)')
		return l
	ret=os.popen(context.env['CXX']+' -print-file-name=libstdc++.so.6').readlines()[0][:-1]
	if ret[0]!='/':
		context.Result('Relative path "%s" was given by compiler %s, must specify libstdcxx=.. explicitly.'%(ret,context.env['CXX']))
		Exit(1)
	ret=os.path.abspath(ret) # removes .. in the path returned by g++
	context.env['libstdcxx']=ret
	context.Result(ret)
	return ret

def CheckPython(context):
	"Checks for functional python/c API. Sets variables if OK and returns true; otherwise returns false."
	origs={'LIBS':context.env['LIBS'],'LIBPATH':context.env['LIBPATH'],'CPPPATH':context.env['CPPPATH'],'LINKFLAGS':context.env['LINKFLAGS']}
	context.Message('Checking for Python development files... ')
	try:
		#FIXME: once caught, exception disappears along with the actual message of what happened...
		import distutils.sysconfig as ds
		context.env.Append(CPPPATH=ds.get_python_inc(),LIBS=ds.get_config_var('LIBS').split() if ds.get_config_var('LIBS') else None)
		context.env.Append(LINKFLAGS=ds.get_config_var('LINKFORSHARED').split()+ds.get_config_var('BLDLIBRARY').split())
		ret=context.TryLink('#include<Python.h>\nint main(int argc, char **argv){Py_Initialize(); Py_Finalize();}\n','.cpp')
		if not ret: raise RuntimeError
	except (ImportError,RuntimeError,ds.DistutilsPlatformError):
		for k in origs.keys(): context.env[k]=origs[k]
		context.Result('error')
		return False
	context.Result('ok'); return True

def CheckBoost(context):
	context.Message('Checking boost libraries... ')
	libs=[
		('boost_system','boost/system/error_code.hpp','boost::system::error_code();',False),
		('boost_thread','boost/thread/thread.hpp','boost::thread();',True),
		('boost_date_time','boost/date_time/posix_time/posix_time.hpp','boost::posix_time::time_duration();',True),
		('boost_filesystem','boost/filesystem/path.hpp','boost::filesystem::path();',True),
		('boost_iostreams','boost/iostreams/device/file.hpp','boost::iostreams::file_sink("");',True),
		('boost_regex','boost/regex.hpp','boost::regex("");',True),
		('boost_serialization','boost/archive/archive_exception.hpp','try{} catch (const boost::archive::archive_exception& e) {};',True),
		('boost_program_options','boost/program_options.hpp','boost::program_options::options_description o;',True),
		('boost_python','boost/python.hpp','boost::python::scope();',True),
	]
	failed=[]
	def checkLib_maybeMT(lib,header,func):
		for LIB in (lib,lib+'-mt'):
			LIBS=env['LIBS'][:]; context.env.Append(LIBS=[LIB])
			if context.TryLink('#include<'+header+'>\nint main(void){'+func+';}','.cpp'): return True
			env['LIBS']=LIBS
		return False
	for lib,header,func,mandatory in libs:
		ok=checkLib_maybeMT(lib,header,func)
		if not ok and mandatory: failed.append(lib)
	if failed: context.Result('Failures: '+', '.join(failed)); return False
	context.Result('all ok'); return True

def CheckPythonModules(context):
	context.Message("Checking for required python modules... ")
	mods=[('IPython','ipython'),('numpy','python-numpy'),('matplotlib','python-matplotlib'),('Xlib','python-xlib')]
	if 'qt4' in context.env['features']: mods.append(('PyQt4.QtGui','python-qt4'))
	failed=[]
	for m,pkg in mods:
		try:
			exec("import %s"%m)
		except ImportError:
			failed.append(m+' (package %s)'%pkg)
	if failed: context.Result('Failures: '+', '.join(failed)); return False
	context.Result('all ok'); return True

	

if not env.GetOption('clean'):
	conf=env.Configure(custom_tests={'CheckLibStdCxx':CheckLibStdCxx,'CheckCXX':CheckCXX,'CheckBoost':CheckBoost,'CheckPython':CheckPython,'CheckPythonModules':CheckPythonModules}, # 'CheckQt':CheckQt
		conf_dir='$buildDir/.sconf_temp',log_file='$buildDir/config.log'
	)
	ok=True
	ok&=conf.CheckCXX()
	if not ok:
			print "\nYour compiler is broken, no point in continuing. See `%s' for what went wrong and use the CXX/CXXFLAGS parameters to change your compiler."%(buildDir+'/config.log')
			Exit(1)
	conf.CheckLibStdCxx()
	ok&=conf.CheckLibWithHeader('pthread','pthread.h','c','pthread_exit(NULL);',autoadd=1)
	ok&=(conf.CheckPython() and conf.CheckCXXHeader(['Python.h','numpy/ndarrayobject.h'],'<>'))
	ok&=conf.CheckPythonModules()
	ok&=conf.CheckBoost()
	env['haveForeach']=conf.CheckCXXHeader('boost/foreach.hpp','<>')
	if not env['haveForeach']: print "(OK, local version will be used instead)"
	ok&=conf.CheckCXXHeader('Eigen/Core')
	ok&=conf.CheckCXXHeader('loki/NullType.h')
	# for installable stript's shebang ( http://en.wikipedia.org/wiki/Shebang_(Unix) )
	env['pyExecutable']=sys.executable

	if not ok:
		print "\nOne of the essential libraries above was not found, unable to continue.\n\nCheck `%s' for possible causes, note that there are options that you may need to customize:\n\n"%(buildDir+'/config.log')+opts.GenerateHelpText(env)
		Exit(1)
	def featureNotOK(featureName,note=None):
		print "\nERROR: Unable to compile with optional feature `%s'.\n\nIf you are sure, remove it from features (scons features=featureOne,featureTwo for example) and build again."%featureName
		if note: print "Note:",note
		Exit(1)
	# check "optional" libs
	if 'opengl' in env['features']:
		ok=conf.CheckLibWithHeader('glut','GL/glut.h','c++','glutGetModifiers();',autoadd=1)
		# TODO ok=True for darwin platform where openGL (and glut) is native
		if not ok: featureNotOK('opengl')
		if 'qt4' in env['features']:
			env['ENV']['PKG_CONFIG_PATH']='/usr/bin/pkg-config'
			env.Tool('qt4')
			env.EnableQt4Modules(['QtGui','QtCore','QtXml','QtOpenGL'])
			if not conf.TryAction(env.Action('pyrcc4'),'','qrc'): featureNotOK('qt4','The pyrcc4 program is not operational (package pyqt4-dev-tools)')
			if not conf.TryAction(env.Action('pyuic4'),'','ui'): featureNotOK('qt4','The pyuic4 program is not operational (package pyqt4-dev-tools)')
			if conf.CheckLibWithHeader(['qglviewer-qt4'],'QGLViewer/qglviewer.h','c++','QGLViewer();',autoadd=1):
				env['QGLVIEWER_LIB']='qglviewer-qt4'
			elif conf.CheckLibWithHeader(['libQGLViewer'],'QGLViewer/qglviewer.h','c++','QGLViewer();',autoadd=1):
				env['QGLVIEWER_LIB']='libQGLViewer'
			else: featureNotOK('qt4','Building with Qt4 implies the QGLViewer library installed (package libqglviewer-qt4-dev package in debian/ubuntu, libQGLViewer in RPM-based distributions)')
	if 'vtk' in env['features']:
		ok=conf.CheckLibWithHeader(['vtkCommon'],'vtkInstantiator.h','c++','vtkInstantiator::New();',autoadd=1)
		env.Append(LIBS=['vtkHybrid','vtkFiltering','vtkRendering','vtkIO','vtkexoIIc','vtkParallel','vtkGraphics','vtkverdict','vtkImaging','vtkftgl','vtkDICOMParser','vtkmetaio','vtksqlite'])
		if not ok: featureNotOK('vtk',note="Installer can`t find vtk-library. Be sure you have it installed (usually, libvtk5-dev package). Or you might have to add VTK header directory (e.g. /usr/include/vtk-5.4) to CPPPATH.")
	if 'gts' in env['features']:
		env.ParseConfig('pkg-config gts --cflags --libs');
		ok=conf.CheckLibWithHeader('gts','gts.h','c++','gts_object_class();',autoadd=1)
		env.Append(CPPDEFINES='PYGTS_HAS_NUMPY')
		if not ok: featureNotOK('gts')
	if 'log4cxx' in env['features']:
		ok=conf.CheckLibWithHeader('log4cxx','log4cxx/logger.h','c++','log4cxx::Logger::getLogger("");',autoadd=1)
		if not ok: featureNotOK('log4cxx')
	if 'gl2ps' in env['features']:
		ok=conf.CheckLibWithHeader('gl2ps','gl2ps.h','c','gl2psEndPage();',autoadd=1)
		if not ok: featureNotOK('gl2ps')
	if 'cgal' in env['features']:
		ok=conf.CheckLibWithHeader('CGAL','CGAL/Exact_predicates_inexact_constructions_kernel.h','c++','CGAL::Exact_predicates_inexact_constructions_kernel::Point_3();')
		env.Append(CXXFLAGS='-frounding-math') # required by cgal, otherwise we get assertion failure at startup
		if not ok: featureNotOK('cgal')
	if not env['mono']: env.Append(LIBS='yade-support')

	env.Append(CPPDEFINES=['YADE_'+f.upper().replace('-','_') for f in env['features']])

	env=conf.Finish()

	if os.path.exists('../brefcom-mm.hh'):
		print "Will use full CPM model in ../brefcom-mm.hh"
		env.Append(CPPDEFINES='YADE_CPM_FULL_MODEL_AVAILABLE')

##########################################################################################
############# BUILDING ###################################################################
##########################################################################################

### SCONS OPTIMIZATIONS
SCons.Defaults.DefaultEnvironment(tools = [])
env.Decider('MD5-timestamp')
env.SetOption('max_drift',5) # cache md5sums of files older than 5 seconds
SetOption('implicit_cache',0) # cache #include files etc.
env.SourceCode(".",None) # skip dotted directories
SetOption('num_jobs',env['jobs'])

### SHOWING OUTPUT
if env['brief']:
	## http://www.scons.org/wiki/HidingCommandLinesInOutput
	env.Replace(CXXCOMSTR='C ${SOURCES}', # → ${TARGET.file}')
		CCOMSTR='C ${SOURCES}',
		SHCXXCOMSTR='C ${SOURCES}', 
		SHCCCOMSTR='C ${SOURCES}', 
		SHLINKCOMSTR='L ${TARGET.file}',
		LINKCOMSTR='L ${TARGET.file}',
		INSTALLSTR='> $TARGET',
		QT_UICCOMSTR='U ${SOURCES}',
		QT_MOCCOMSTR='M ${SOURCES}',
		QT4_UICCOMSTR='U ${SOURCES}',
		QT_MOCFROMHCOMSTR='M ${SOURCES}',
		QT_MOCFROMCXXCOMSTR='M ${SOURCES}',
		)
else:
	env.Replace(INSTALLSTR='cp -f ${SOURCE} ${TARGET}')

### DIRECTORIES
## PREFIX must be absolute path. Why?!
env['PREFIX']=os.path.abspath(env['PREFIX'])

libDirs=('extra','gui','lib','py','plugins')
# where are we going to be installed... pkg/dem becomes pkg-dem
instLibDirs=[os.path.join('$LIBDIR',x) for x in libDirs]
## runtime library search directories; there can be up to 2 levels of libs, so we do in in quite a crude way here:
## FIXME: use syntax as shown here: http://www.scons.org/wiki/UsingOrigin
relLibDirs=['../'+x for x in libDirs]+['../../'+x for x in libDirs]+[env.subst('../lib/yade$SUFFIX_DBG/lib')]
runtimeLibDirs=[env.Literal('\\$$ORIGIN/'+x) for x in relLibDirs]

### PREPROCESSOR FLAGS
if env['QUAD_PRECISION']: env.Append(CPPDEFINES='QUAD_PRECISION')

### COMPILER
if env['debug']: env.Append(CXXFLAGS='-ggdb2',CPPDEFINES=['YADE_DEBUG'])
if 'openmp' in env['features']: env.Append(CXXFLAGS='-fopenmp',LIBS='gomp',CPPDEFINES='YADE_OPENMP')
if env['optimize']:
	# NDEBUG is used in /usr/include/assert.h: when defined, asserts() are no-ops
	env.Append(CXXFLAGS=['-O3'],CPPDEFINES=[('YADE_CAST','static_cast'),('YADE_PTR_CAST','static_pointer_cast'),'NDEBUG'])
	# do not state architecture if not provided
	# used for debian packages, where 'native' would generate instructions outside the package architecture
	# (such as SSE instructions on i386 builds, if the builder supports them)
	if env.has_key('march') and env['march']: env.Append(CXXFLAGS=['-march=%s'%env['march']]),
else:
	env.Append(CPPDEFINES=[('YADE_CAST','dynamic_cast'),('YADE_PTR_CAST','dynamic_pointer_cast')])

if env['gprof']: env.Append(CXXFLAGS=['-pg'],LINKFLAGS=['-pg'],SHLINKFLAGS=['-pg'])
env.Prepend(CXXFLAGS=['-pipe','-Wall'])

if env['PGO']=='gen': env.Append(CXXFLAGS=['-fprofile-generate'],LINKFLAGS=['-fprofile-generate'])
if env['PGO']=='use': env.Append(CXXFLAGS=['-fprofile-use'],LINKFLAGS=['-fprofile-use'])

if 'clang' in env['CXX']:
	print 'Looks like we use clang, adding some flags to avoid warning flood.'
	env.Append(CXXFLAGS=['-Wno-unused-variable','-Wno-mismatched-tags','-Wno-constant-logical-operand','-Qunused-arguments','-Wno-empty-body'])
	if 'openmp' in env['features']: print 'WARNING: building with clang and OpenMP feature, expect compilation errors!'
	if env['march']: print 'WARNING: specifying march with clang\'s might lead to crash at startup (if so, recompile with march= )!'


### LINKER
## libs for all plugins
if not env['mono']:
	env.Append(LIBS=[],SHLINKFLAGS=['-rdynamic'])
	env.Append(LINKFLAGS=['-rdynamic','-Wl,-z,origin'])

if not env['debug']: env.Append(SHLINKFLAGS=['-W,--strip-all'])

# makes dynamic library loading easier (no LD_LIBRARY_PATH) and perhaps faster
env.Append(RPATH=runtimeLibDirs)
# find already compiled but not yet installed libraries for linking
env.Append(LIBPATH=instLibDirs) # this is if we link to libs that are installed, which is the case now


def relpath(pf,pt):
	"""Returns relative path from pf (path from) to pt (path to) as string.
	Last component of pf MUST be filename, not directory name. It can be empty, though. Legal pf's: 'dir1/dir2/something.cc', 'dir1/dir2/'. Illegal: 'dir1/dir2'."""
	from os.path import sep,join,abspath,split
	apfl=abspath(split(pf)[0]).split(sep); aptl=abspath(pt).split(sep); i=0
	while apfl[i]==aptl[i] and i<min(len(apfl),len(aptl))-1: i+=1
	return sep.join(['..' for j in range(0,len(apfl)-i)]+aptl[i:])

# 1. symlink all headers to buildDir/include before the actual build
# 2. (unused now) instruct scons to install (not symlink) all headers to the right place as well
# 3. set the "install" target as default (if scons is called without any arguments), which triggers build in turn
# Should be cleaned up.

if not env.GetOption('clean'):
	# how to make that executed automatically??! For now, run always.
	#env.AddPreAction(installAlias,installHeaders)
	from os.path import join,split,isabs,isdir,exists,lexists,islink,isfile,sep
	#installHeaders() # install to buildDir always
	#if 0: # do not install headers
	#	installHeaders(env.subst('$PREFIX')) # install to $PREFIX if specifically requested: like "scons /usr/local/include"

	yadeInc=join(buildDir,'include/yade')
	if not os.path.exists(yadeInc): os.makedirs(yadeInc)
	import glob
	for p in ['core','extra']+glob.glob('lib/*')+glob.glob('pkg/*')+glob.glob('gui'):
		link=yadeInc+'/'+p.replace('/','-')
		if os.path.isdir(p) and not os.path.exists(link):
			if lexists(link): os.remove(link) # dangling symlink
			os.symlink(relpath(link,p),link)
	boostDir=buildDir+'/include/boost'
	if not exists(boostDir): os.makedirs(boostDir)
	def mkSymlink(link,target):
		if not exists(link):
			if lexists(link): os.remove(link) # remove dangling symlink
			os.symlink(relpath(link,target),link)
	if not env['haveForeach']:
		mkSymlink(boostDir+'/foreach.hpp','extra/foreach.hpp_local')
	#mkSymlink(boostDir+'/python','py/3rd-party/boost-python-indexing-suite-v2-noSymlinkHeaders')
	mkSymlink(buildDir+'/include/indexing_suite','py/3rd-party/boost-python-indexing-suite-v2-noSymlinkHeaders')
	mkSymlink(boostDir+'/math','extra/floating_point_utilities_v3/boost/math')
	env.Default(env.Alias('install',['$PREFIX/bin','$PREFIX/lib'])) # build and install everything that should go to instDirs, which are $PREFIX/{bin,lib} (uses scons' Install)

env.Export('env');

######
### combining builder
#####
# since paths are aboslute, sometime they are not picked up
# hack it away by writing the combined files into a text file
# and when it changes, force rebuild of all combined files
# changes to $buildDir/combined-files are cheked at the end of 
# the SConscript file
import md5
combinedFiles=env.subst('$buildDir/combined-files')

if os.path.exists(combinedFiles):
	combinedFilesMd5=md5.md5(file(combinedFiles).read()).hexdigest()
	os.remove(combinedFiles)
else:
	combinedFilesMd5=None

# http://www.scons.org/wiki/CombineBuilder
import SCons.Action
import SCons.Builder
def combiner_build(target, source, env):
	"""generate a file, that's including all sources"""
	out=""
	for src in source: out+="#include \"%s\"\n"%src.abspath
	open(str(target[0]),'w').write(out)
	env.AlwaysBuild(target[0])
	return 0
def CombineWrapper(target,source):
	open(combinedFiles,'a').write(env.subst(target)+': '+' '.join([env.subst(s) for s in source])+'\n')
	return env.Combine(target,source)
env.CombineWrapper=CombineWrapper
	
env.Append(BUILDERS = {'Combine': env.Builder(action = SCons.Action.Action(combiner_build, "> $TARGET"),target_factory = env.fs.File,)})

import yadeSCons
allPlugs=yadeSCons.scanAllPlugins(None,feats=env['features'])
buildPlugs=yadeSCons.getWantedPlugins(allPlugs,env['exclude'],env['features'],env['chunkSize'],env['hotPlugins'].split(','))
def linkPlugins(plugins):
	"""Given list of plugins we need to link to, return list of real libraries that we should link to."""
	ret=set()
	for p in plugins:
		if not buildPlugs.has_key(p):
			raise RuntimeError("Plugin %s is required (backtrace shows where), but will not be built!"%p)
		ret.add(buildPlugs[p].obj)
	return ['core','yade-support']+list(ret)

env['linkPlugins']=linkPlugins
env['buildPlugs']=buildPlugs

if env['mono']:
	env['topLevelDir']=os.path.abspath(os.getcwd())
	env.SConscript('SConscript-mono',variant_dir=buildDir,duplicate=0)

else:
	# read top-level SConscript file. It is used only so that variant_dir is set. This file reads all necessary SConscripts
	env.SConscript(dirs=['.'],variant_dir=buildDir,duplicate=0)

#################################################################################
## remove plugins that are in the target dir but will not be installed now
## only when installing without requesting special path (we would have no way
## to know what should be installed overall.
if not COMMAND_LINE_TARGETS:
	toInstall=set([str(node) for node in env.FindInstalledFiles()])
	for root,dirs,files in os.walk(env.subst('$LIBDIR')):
		# do not go inside the debug directly, plugins are different there
		for f in files:
			# skip debug files, if in the non-debug build
			if not env['debug'] and '/dbg/' in root: continue
			#print 'Considering',f
			ff=os.path.join(root,f)
			# do not delete python-optimized files and symbolic links (lib_gts__python-module.so, for instance)
			if ff not in toInstall and not ff.endswith('.pyo') and not ff.endswith('.pyc') and not os.path.islink(ff) and not os.path.basename(ff).startswith('.nfs'):
				print "Deleting extra plugin", ff
				os.remove(ff)

#################################################################################
#### DOCUMENTATION
# must be explicitly requested to be installed, e.g.:
#    scons /usr/local/share/doc
env.Install('$PREFIX/share/doc/yade$SUFFIX-doc/',['examples','scripts','doc'])


### check if combinedFiles is different; if so, force rebuild of all of them
if os.path.exists(combinedFiles) and combinedFilesMd5!=md5.md5(open(combinedFiles).read()).hexdigest():
	print 'Rebuilding combined files, since the md5 has changed.'
	combs=[l.split(':')[0] for l in open(combinedFiles)]
	for c in combs:
		env.AlwaysBuild(c)
		env.Default(c)

#Progress('.', interval=100, file=sys.stderr)
