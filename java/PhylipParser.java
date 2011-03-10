/*
	A class to read in sequence data adhering to the Phylip format, as defined
	at http://www.ncbi.nlm.nih.gov/blast/fasta.shtml
*/

import java.util.*;
import java.io.*;
import java.lang.*;

public class PhylipParser
{
	LinkedList<Sequence> sequences;
	BufferedReader input;
	char form; // whether the file is sequential or interleaved.
	
	public PhylipParser(File file, char f)
	{
		sequences = new LinkedList<Sequence>();
		form = f; 
		
		try {
			input = new BufferedReader(new FileReader(file));
		}
		catch(FileNotFoundException fnfe) {
			System.err.println("file " + file.getName() + " not located, exiting");
			fnfe.printStackTrace();
			System.exit(1);
		}
	}
	public PhylipParser(String filename, char f) {
		sequences = new LinkedList<Sequence>();
		form = f;
		
		//instantiates the BufferedReader
		try {
			input = new BufferedReader(new FileReader(new File(filename)));
		}
		catch(FileNotFoundException fnfe) {
			System.err.println("file "+filename+" not located, exiting");
			fnfe.printStackTrace();
			System.exit(1);
		}
	}
	
	/*
		Parses the file, collecting and storing the sequence and description in a new Sequence.
	*/
	public LinkedList<Sequence> parseFile() 
	{
		String currstr = null;
		
		// read in first line of file
		try {
			currstr = input.readLine();
		}
		catch(IOException ioe) {
			System.err.println("ERROR encountered reading in line, exiting");
			ioe.printStackTrace();
			System.exit(1);
		}
		
		// First line of a Phylip file contains the number of taxa and the length of each sequence
		StringTokenizer st = new StringTokenizer(currstr);
		int taxa  = new Integer(st.nextToken());
		int seql = new Integer(st.nextToken());

			
		try {
			currstr = (input.readLine()).trim();
		}
		catch(IOException ioe) {
			System.err.println("ERROR encountered reading in line, exiting");
			ioe.printStackTrace();
			System.exit(1);
		}

		// If the file is in the interleaved format
		if(form == 'I')
		{
			int counttaxa = 0;
			boolean finished = false;
			String tempstr = new String();
			String descrip = new String();
			ArrayList<String> lines = new ArrayList<String>(); // will store all the lines from the file
			
			// Step 1: Read in all the non-blank lines from the file
			while(currstr != null)
			{
				currstr = currstr.trim();
				if(currstr.length() > 1){
					currstr = currstr.toUpperCase();
					lines.add(currstr);
				}
				try {
					currstr = (input.readLine());
				}
				catch(IOException ioe) {
					System.err.println("ERROR encountered reading in line, exiting");
					ioe.printStackTrace();
					System.exit(1);
				}
	
			}
			
			// Step 2: Determine how many lines each sequence covers
			int numlines = lines.size()/taxa;
		
			// Step 3: For each taxa, concat each of the lines together to form the sequence.
			// For example, if taxa1 had the first sequence in the file and there where 10 taxa in the file, concat lines 1,11,21...numlines*10+1
			while(counttaxa < taxa)
			{
				for(int i = 0; i< numlines; i++)
				{
					String line = lines.get(taxa*i + counttaxa);
					String tempshort = new String();
					if(i==0){
						int start = line.indexOf("	");
						if(start >= line.length() || start < 2){
							start = line.indexOf(" ");
						}
						descrip = line.substring(0,start);
						tempshort = removeSpace(line.substring(start+1));
					}
					else{
						tempshort = removeSpace(line);
					}
					tempstr += tempshort;
				}
		
				Sequence s = new Sequence(tempstr,descrip);
				tempstr = "";
				sequences.add(s);
				counttaxa++;
			}
		}
		// If the file is in the sequential format
		else if(form == 'S')
		{
			String tempstr = new String();
			String descrip = new String();
			
			// Step 1: Determine how many bases/aminoacids there are per line and hence how many lines per sequence

			int start = currstr.indexOf("	");
			if(start >= currstr.length() || start < 2){
				start = currstr.indexOf(" ");
			}

			String firstline  = removeSpace(currstr.substring(start+1).trim());
			int bpl = firstline.length(); // no of bases/amino acids per line;
			int lines=0;
			if(seql==bpl){lines=1;}
			else{
			lines = seql/bpl + 1;}
			
			// Step 2: For each taxa, concat each of the lines together to form the sequence.
			int count = 1;
			while(currstr!=null) 
			{
				if(currstr.length() > 0){
					currstr = (currstr.toUpperCase()).trim();
				}
				if(count == lines && lines == 1)
				{
					start = currstr.indexOf(" ");
					if(start >= currstr.length() || start < 2){
						start = currstr.indexOf("	");
					}
					descrip = currstr.substring(0,start);
					tempstr = removeSpace(currstr.substring(start+1));
					Sequence seq = new Sequence(tempstr,descrip);
					sequences.addLast(seq);
					count = 0;
				}
				else if(count == lines && lines > 1)
				{
					tempstr+=removeSpace(currstr);
					Sequence seq = new Sequence(tempstr,descrip);
					sequences.addLast(seq);
					count = 0;
				}
				else if(count == 1)
				{
					start = currstr.indexOf(" ");
					if(start >= currstr.length() || start < 2){
						start = currstr.indexOf("	");
					}
					descrip = currstr.substring(0,start);
					tempstr = removeSpace(currstr.substring(start+1));
					tempstr = tempstr;
				}
				else
				{
					tempstr += removeSpace(currstr);
				}
				count++;
				
				//moves onto next line of input file
				try {
					currstr = input.readLine();
				}
				catch(IOException ioe) {
					System.err.println("ERROR encountered reading in line, exiting");
					ioe.printStackTrace();
					System.exit(1);
				}
			}
		}
		
		return sequences;
	}
	
	/*
		Method to remove any white space contained within the sequence.
		Eg. If characters in sequence are grouped in tens, with a space separating the groups, this method will remove the spaces.
	*/
	public String removeSpace(String oldString)
	{
		String str = new String();
		String[] newString = oldString.split(" ");
		for(int i = 0; i <  newString.length; i++)
		{
			str += newString[i];
		}
		return str;
	}
}