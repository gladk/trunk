#!/usr/bin/python
# encoding: utf-8
try:
	import _bibtex as bib
	import _recode as bibRecode
except ImportError:
	raise ImportError("Unable to import _bibtex and/or _recode; install the python-bibtex package.")
import sys

def readBib(filename):
	## _bibtex has horrible interface
	bibfile=bib.open_file(filename,1) # 2nd argument: strict mode
	db={}
	#rq=bibRecode.request('latex..latin1')
	while True:
		entry=bib.next(bibfile)
		if entry is None: break
		key,type,dta=entry[0],entry[1],entry[4]
		item={'type':type}
		for field in dta.keys():
			expanded=bib.expand(bibfile,dta[field],-1)
			#conv=bibRecode.recode(rq,expanded[2])
			item[field]=expanded[2].strip()
		db[key]=item
	## now we don't need _bibtex anymore, everything is in our dicts
	return db

def dumpBib(db):
	for k in db.keys():
		print k,db[k]

def formatRest(db):
	ret=[]
	keys=db.keys(); keys.sort()
	for key in keys:
		i=db[key]; type=i['type']
		line=r'.. [%s] \ '%key ## ← HACK: explicit space to prevent docutils from using abbreviated first name (e.g. "F.") as enumeration item; it works!!
		#ret.append('.. %s:'%key)
		#ret.append(':[%s]:'%key)
		#line='\t\ ' ## this also prevent considering abbreviated initial as enumeration

		# required fields from http://en.wikipedia.org/wiki/Bibtex
		# in some cases, they are not present, anyway
		if type=='article':
			if i.has_key('author'): line+='%s '%i['author']
			if i.has_key('year'): line+='(%s), '%i['year']
			line+='**%s**. *%s*'%(i['title'],i['journal'])
			if i.has_key('issue'): line+=i['issue']
			if i.has_key('volume'): line+=' (%s)'%i['volume']
			if i.has_key('pages'): line+=', pages %s'%i['pages']
			line+='.'
		elif type=='book':
			if i.has_key('author'): line+='%s '%i['author']
			if i.has_key('year'): line+='(%s), '%i['year']
			line+='**%s**.'%i['title']
			if i.has_key('publisher'): line+=' %s.'%i['publisher']
		elif type=='inproceedings':
			line+='%s (%s), **%s**. In *%s*.'%(i['author'],i['year'],i['title'],i['booktitle'])
		elif type=='phdthesis':
			line+='%s (%s), **%s**. PhD thesis at *%s*.'%(i['author'],i['year'],i['title'],i['school'])
		elif type=='proceedings':
			if i.has_key('editor'): line+='%s (ed.), '%i['editor']
			line+='**%s** (%s).'%(i['title'],i['year'])
			if i.has_key('organization'): line+=' *%s*.'%i['organization']
			if i.has_key('publisher'): line+=' %s'%i['publisher']
		elif type=='misc':
			if i.has_key('author'): line+=i['author']
			if i.has_key('year'): line+=' %s'%i['year']
			if i.has_key('title'): line+=' **%s**'%i['title']
		# add doi and url to everything, if present
		## ReST uses <..> to delimit URL, therefore < and > must be encoded in the URL (http://www.blooberry.com/indexdot/html/topics/urlencoding.htm)
		if i.has_key('doi'): line+=' DOI `%s <http://dx.doi.org/%s>`_'%(i['doi'],i['doi'].replace('<','%3c').replace('>','%3e')) 
		if i.has_key('url'): line+=' URL %s'%i['url']
		ret.append(line)
	return ret

def bib2rst(filename):
	"""Return string representing all items in given bibtex file, formatted as ReST."""
	db=readBib(filename)
	return '\n\n'.join(formatRest(db))
if __name__=='__main__':
	import sys
	print bib2rst(sys.argv[1])
