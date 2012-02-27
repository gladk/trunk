# encoding: utf-8
# 2009 © Václav Šmilauer <eudoxos@arcig.cz>
"""All defined functionality tests for yade."""
import unittest,inspect

# add any new test suites to the list here, so that they are picked up by testAll
allTests=['wrapper','omega','pbc','clump'] # ,'cohesive-chain']

# all yade modules (ugly...)
import yade.eudoxos,yade.export,yade.linterpolation,yade.log,yade.pack,yade.plot,yade.post2d,yade.timing,yade.utils,yade.ymport
allModules=(yade.eudoxos,yade.export,yade.linterpolation,yade.log,yade.pack,yade.plot,yade.post2d,yade.timing,yade.utils,yade.ymport)
try:
	import yade.qt
	allModules+=(yade.qt,)
except ImportError: pass

# fully qualified module names
allTestsFQ=['yade.tests.'+test for test in allTests]

def testModule(module):
	"""Run all tests defined in the module specified, return TestResult object 
	(http://www.python.org/doc/2.6/lib/unittest.html#testresult-objects.html)
	for further processing.

	@param module: fully-qualified module name, e.g. yade.tests.wrapper
	"""
	suite=unittest.defaultTestLoader().loadTestsFromName(module)
	return unittest.TextTestRunner(verbosity=2).run(suite)

def testAll():
	"""Run all tests defined in all yade.tests.* modules and return
	TestResult object for further examination."""
	suite=unittest.defaultTestLoader.loadTestsFromNames(allTestsFQ)
	import doctest
	for mod in allModules:
		suite.addTest(doctest.DocTestSuite(mod))
	return unittest.TextTestRunner(verbosity=2).run(suite)

	

