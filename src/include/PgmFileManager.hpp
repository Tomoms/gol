#ifndef PGMFILEMANAGER_H
#define PGMFILEMANAGER_H

#include <string>

#define PGM_MAX_VALUE	255
#define PGM_HOLDER		std::vector<std::vector<unsigned char>>

class PgmFileManager
{

private:
	std::string filename_;
	unsigned long rows_, cols_;

public:
	PgmFileManager(std::string& filename); // open existing file ctor
	PgmFileManager(std::string& filename, unsigned long size); // new file ctor
	void write();
	PGM_HOLDER read();
	unsigned long get_rows();
	unsigned long get_cols();
};

#endif
