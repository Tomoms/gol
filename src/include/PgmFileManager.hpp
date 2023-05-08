#ifndef PGMFILEMANAGER_H
#define PGMFILEMANAGER_H

#define PGM_MAX_VALUE	255

class PgmFileManager
{

private:
	std::string filename;
	unsigned long rows, cols;

public:
	//PgmFileManager(std::string& filename); // open existing file ctor
	PgmFileManager(std::string& filename, unsigned long size); // new file ctor
	void write();
};

#endif
