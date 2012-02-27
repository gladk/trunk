#pragma once
#include<yade/core/FileGenerator.hpp>
#include<yade/pkg-dem/Shop.hpp>


class SimpleScene: public FileGenerator {
	public:
		SimpleScene(){};
		~SimpleScene (){};
		virtual bool generate();
	REGISTER_ATTRIBUTES(FileGenerator,/* */);
	REGISTER_CLASS_NAME(SimpleScene);
	REGISTER_BASE_CLASS_NAME(FileGenerator);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(SimpleScene);

