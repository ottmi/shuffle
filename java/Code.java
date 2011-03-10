import java.util.*;
import java.lang.*;

/*
	Class for storage of information on how program interprets the sequence data supplied to it
	Handles the type of code (nucleotide, amino acid, user defined taxa) and sequence character
	groupings (one character at a time, two, three ...).
*/
public class Code
{
	private String name;
	private LinkedList<String> members;	//linkedlist containing all possible discrete sequence characters
	private char unknown;	//the character that represents an unknown sequence character
	private char missing;
	private int character_grouping_number;	//the number of discrete characters in each site member
	private HashMap<String,String> nucleotide_labels;
	private HashMap<String,String[]> amino_acid_labels;
	private HashMap<String,String> codon_translations;
	private int sequence_data_type;
	private int offset_number;	//the number of discrete characters, from the beginning of a sequence, the program should skip before starting to create sites
	private HashSet<String> character_groups;
	
	
	//static variables defining the type of sequence data - used at construction
	public static final int AMINO_ACID = 0;
	public static final int NUCLEOTIDE_DNA = 1;
	public static final int NUCLEOTIDE_RNA = 2;
	public static final int CUSTOM_TAXA = 3;
	public static final int STANDARD = 4;
	
	
	public static final String NUCLEOTIDE_CHARACTERS_CERTAIN = "ATGCU";
	public static final String NUCLEOTIDE_CHARACTERS_AMBIGUOUS = "MRWSYKVHDBN";
	public static final String NUCLEOTIDE_CHARACTERS_ALL = "ATGCUMRWSYKVHDBN";
	
	public static final String AMINO_ACID_CHARACTERS = "ABCDEFGHIKLMNPQRSTUVWYZX*";
	
	public static final String STANDARD_CHARACTERS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789~!@#$%^&*+=";
	
	public static final String[] AMINO_ACID_CHARACTERS_SINGLE = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N",
																	"P", "Q", "R", "S", "T", "U", "V", "W", "Y"};
	public static final String[] STANDARD_CHARACTERS_SINGLE = {"A", "B","C", "D", "E", "F", "G", "H", "I", "J","K","L", "M", "N","O",
																	"P", "Q", "R", "S", "T", "U", "V", "W", "X","Y","Z","0","1","2","3","4","5","6","7","8","9",
																			"~","!","@","#","$","%","^","&","*","+","="};
	
	public static final String[] NUCLEOTIDE_CHARACTERS_SINGLE_DNA = {"A", "T", "G", "C"};
	
	public static final String[] NUCLEOTIDE_CHARACTERS_SINGLE_RNA = {"A", "U", "G", "C"};
	
	public static final String[] NUCLEOTIDE_CHARACTERS_SINGLE_AMBIGUOUS = {"M","R","W","S","Y","K","V","H","D","B","N"};
	
	public static final String[] NUCLEOTIDE_CHARACTERS_DOUBLE_DNA = {"AA", "AT", "AG", "AC", 
																		"TA", "TT", "TG", "TC",
																		"GA", "GT", "GG", "GC",
																		"CA", "CT", "CG", "CC"};
	public static final String[] NUCLEOTIDE_CHARACTERS_DOUBLE_RNA = {"AA", "AU", "AG", "AC", 
																		"UA", "UU", "UG", "UC",
																		"GA", "GU", "GG", "GC",
																		"CA", "CU", "CG", "CC"};
	
	public static final String[] NUCLEOTIDE_CHARACTERS_TRIPLE_DNA = {"AAA", "AAT", "AAG", "AAC", 
																		"ATA", "ATT", "ATG", "ATC",
																		"AGA", "AGT", "AGG", "AGC",
																		"ACA", "ACT", "ACG", "ACC",																		
																		"TAA", "TAT", "TAG", "TAC", 
																		"TTA", "TTT", "TTG", "TTC",
																		"TGA", "TGT", "TGG", "TGC",
																		"TCA", "TCT", "TCG", "TCC",																		
																		"GAA", "GAT", "GAG", "GAC", 
																		"GTA", "GTT", "GTG", "GTC",
																		"GGA", "GGT", "GGG", "GGC",
																		"GCA", "GCT", "GCG", "GCC",																		
																		"CAA", "CAT", "CAG", "CAC", 
																		"CTA", "CTT", "CTG", "CTC",
																		"CGA", "CGT", "CGG", "CGC",
																		"CCA", "CCT", "CCG", "CCC"};
	
	public static final String[] NUCLEOTIDE_CHARACTERS_TRIPLE_RNA = {"AAA", "AAU", "AAG", "AAC", 
																		"AUA", "AUU", "AUG", "AUC",
																		"AGA", "AGU", "AGG", "AGC",
																		"ACA", "ACU", "ACG", "ACC",																		
																		"UAA", "UAU", "UAG", "UAC", 
																		"UUA", "UUU", "UUG", "UUC",
																		"UGA", "UGU", "UGG", "UGC",
																		"UCA", "UCU", "UCG", "UCC",																		
																		"GAA", "GAU", "GAG", "GAC", 
																		"GUA", "GUU", "GUG", "GUC",
																		"GGA", "GGU", "GGG", "GGC",
																		"GCA", "GCU", "GCG", "GCC",																		
																		"CAA", "CAU", "CAG", "CAC", 
																		"CUA", "CUU", "CUG", "CUC",
																		"CGA", "CGU", "CGG", "CGC",
																		"CCA", "CCU", "CCG", "CCC"};


	
	public String nucleotide_characters;
	public String amino_acid_characters;
	
	public Code(String codeName, int sequenceType, char unknownCharacter, char missingCharacter, int groupingNumber)
	{
		this(codeName, sequenceType, unknownCharacter,missingCharacter, groupingNumber, 0);
		
	}
	public Code(String codeName, int sequenceType, char unknownCharacter, char missingCharacter, int groupingNumber, int offsetNumber) {
		
		name = codeName;
		members = new LinkedList<String>();
		unknown = unknownCharacter;
		missing = missingCharacter;
		character_grouping_number = groupingNumber;
		sequence_data_type = sequenceType;
		offset_number = offsetNumber;
		
		initialiseHashMaps();
		setCharacterFields();
	}
	
	
	/*
		sets the character fields to their appropriate values, using those
		defined in the hashmaps appropriate to the type. must always be called
		at construction AFTER initialiseHashMaps()
	*/
	private void setCharacterFields() {
		//set nucleotides
		Set<String> nuc_set = nucleotide_labels.keySet();
		String[] nuc_labels = nuc_set.toArray(new String[1]);
		nucleotide_characters = new String();
		for(int i = 0; i < nuc_labels.length; i++) { nucleotide_characters = nucleotide_characters.concat(nuc_labels[i]); }
		
		//set amino acids
		Set<String> aa_set = amino_acid_labels.keySet();
		String[] aa_labels = aa_set.toArray(new String[1]);
		amino_acid_characters = new String();
		for(int i = 0; i < aa_labels.length; i++) { amino_acid_characters = amino_acid_characters.concat(aa_labels[i]); }
	}
	
	/* 
		puts all the mappings into appropriate HashMaps. at construction, must always be called
		before any other method relying on instantiated hashmaps
	*/
	private void initialiseHashMaps() {
		nucleotide_labels = new HashMap<String,String>(16);
		nucleotide_labels.put("A", "Adenine");
		nucleotide_labels.put("T", "Thymine");
		nucleotide_labels.put("G", "Guanine");
		nucleotide_labels.put("C", "Cytosine");
		nucleotide_labels.put("U", "Uracil");
		nucleotide_labels.put("R", "Guanine/Adenine (Purine)");
		nucleotide_labels.put("Y", "Thymine/Cytosine (Pyrimidine)");
		nucleotide_labels.put("K", "Guanine/Thymine (Ketone)");
		nucleotide_labels.put("M", "Adenine/Cytosine (Amino group)");
		nucleotide_labels.put("S", "Guanine/Cytosine (Strong interaction)");
		nucleotide_labels.put("W", "Adenine/Thymine (Weak interaction)");
		nucleotide_labels.put("B", "Guanine/Thymine/Cytosine (Not Adenine)");
		nucleotide_labels.put("D", "Guanine/Adenine/Thymine (Not Cytosine)");
		nucleotide_labels.put("H", "Adenine/Cytosine/Thymine (Not Guanine)");
		nucleotide_labels.put("V", "Guanine/Cytosine/Adenine (Not Thymine/Uracil)");
		nucleotide_labels.put("N", "Adenine/Guanine/Cytosine/Thymine");
	
		
		amino_acid_labels = new HashMap<String,String[]>(25);
		amino_acid_labels.put("A", new String[] {"Ala", "Alanine"});
		amino_acid_labels.put("B", new String[] {"Asp/Asn", "Aspartic Acid/Asparagine"});
		amino_acid_labels.put("C", new String[] {"Cys", "Cysteine"});
		amino_acid_labels.put("D", new String[] {"Asp", "Aspartic Acid"});
		amino_acid_labels.put("E", new String[] {"Glu", "Glutamic Acid"});
		amino_acid_labels.put("F", new String[] {"Phe", "Phenylalanine"});
		amino_acid_labels.put("G", new String[] {"Gly", "Glycine"});
		amino_acid_labels.put("H", new String[] {"His", "Histidine"});
		amino_acid_labels.put("I", new String[] {"Ile", "Isoleucine"});
		amino_acid_labels.put("K", new String[] {"Lys", "Lysine"});
		amino_acid_labels.put("L", new String[] {"Leu", "Leucine"});
		amino_acid_labels.put("M", new String[] {"Met", "Methionine"});
		amino_acid_labels.put("N", new String[] {"Asn", "Asparagine"});
		amino_acid_labels.put("P", new String[] {"Pro", "Proline"});
		amino_acid_labels.put("Q", new String[] {"Gln", "Glutamine"});
		amino_acid_labels.put("R", new String[] {"Arg", "Arginine"});
		amino_acid_labels.put("S", new String[] {"Ser", "Serine"});
		amino_acid_labels.put("T", new String[] {"Thr", "Threonine"});
		amino_acid_labels.put("U", new String[] {"Sec", "Selenocysteine"});
		amino_acid_labels.put("V", new String[] {"Val", "Valine"});
		amino_acid_labels.put("W", new String[] {"Trp", "Tryptophan"});
		amino_acid_labels.put("Y", new String[] {"Tyr", "Tyrosine"});
		amino_acid_labels.put("Z", new String[] {"Glu/Gln", "Glutamic Acid/Glutamine"});
		amino_acid_labels.put("X", new String[] {"Any", "Any"});
		amino_acid_labels.put("*", new String[] {"End", "Stop/Termination"});
		
		
		codon_translations = new HashMap<String,String>(64);
		codon_translations.put("TTT", "F");
		codon_translations.put("TTC", "F");
		codon_translations.put("TTA", "L");
		codon_translations.put("TTG", "L");
		codon_translations.put("TCT", "S");
		codon_translations.put("TCC", "S");
		codon_translations.put("TCA", "S");
		codon_translations.put("TCG", "S");
		codon_translations.put("TAT", "Y");
		codon_translations.put("TAC", "Y");
		codon_translations.put("TAA", "*");
		codon_translations.put("TAG", "*");
		codon_translations.put("TGT", "C");
		codon_translations.put("TGC", "C");
		codon_translations.put("TGA", "*");
		codon_translations.put("TGG", "W");
		
		codon_translations.put("CTT", "L");
		codon_translations.put("CTC", "L");
		codon_translations.put("CTA", "L");
		codon_translations.put("CTG", "L");
		codon_translations.put("CCT", "P");
		codon_translations.put("CCC", "P");
		codon_translations.put("CCA", "P");
		codon_translations.put("CCG", "P");
		codon_translations.put("CAT", "H");
		codon_translations.put("CAC", "H");
		codon_translations.put("CAA", "Q");
		codon_translations.put("CAG", "Q");
		codon_translations.put("CGT", "R");
		codon_translations.put("CGC", "R");
		codon_translations.put("CGA", "R");
		codon_translations.put("CGG", "R");
		
		codon_translations.put("ATT", "I");
		codon_translations.put("ATC", "I");
		codon_translations.put("ATA", "I");
		codon_translations.put("ATG", "M");
		codon_translations.put("ACT", "T");
		codon_translations.put("ACC", "T");
		codon_translations.put("ACA", "T");
		codon_translations.put("ACG", "T");
		codon_translations.put("AAT", "N");
		codon_translations.put("AAC", "N");
		codon_translations.put("AAA", "K");
		codon_translations.put("AAG", "K");
		codon_translations.put("AGT", "S");
		codon_translations.put("AGC", "S");
		codon_translations.put("AGA", "R");
		codon_translations.put("AGG", "R");
		
		codon_translations.put("GTT", "V");
		codon_translations.put("GTC", "V");
		codon_translations.put("GTA", "V");
		codon_translations.put("GTG", "V");
		codon_translations.put("GCT", "A");
		codon_translations.put("GCC", "A");
		codon_translations.put("GCA", "A");
		codon_translations.put("GCG", "A");
		codon_translations.put("GAT", "D");
		codon_translations.put("GAC", "D");
		codon_translations.put("GAA", "E");
		codon_translations.put("GAG", "E");
		codon_translations.put("GGT", "G");
		codon_translations.put("GGC", "G");
		codon_translations.put("GGA", "G");
		codon_translations.put("GGG", "G");
	}
	
	
	//get methods
	public String getName()	{ return name; }
	public int getSequenceType() { return sequence_data_type; }
	
	
	//will only work right if determinedContainedMembers() has been previously called
	public String[] getMembers()
	{
		return character_groups.toArray(new String[] {});
	}
	
	
	//returns arrays containing possible site contents
	public String[] getMembers2() {
		if(sequence_data_type == Code.NUCLEOTIDE_DNA || sequence_data_type == Code.NUCLEOTIDE_RNA) 
		{
			if(sequence_data_type == Code.NUCLEOTIDE_DNA) {
				if(character_grouping_number == 1) {
					return Code.NUCLEOTIDE_CHARACTERS_SINGLE_DNA;
				}
				else if(character_grouping_number == 2) {
					return Code.NUCLEOTIDE_CHARACTERS_DOUBLE_DNA;
				}
				else if(character_grouping_number == 3) {
					return Code.NUCLEOTIDE_CHARACTERS_TRIPLE_DNA;
				}
				else {
					System.out.println("ERROR: Too high a grouping number! Exiting");
					System.exit(1);
					return null;
				}
			}
			else if(sequence_data_type == Code.NUCLEOTIDE_RNA) {
				if(character_grouping_number == 1) {
					return Code.NUCLEOTIDE_CHARACTERS_SINGLE_RNA;
				}
				else if(character_grouping_number == 2) {
					return Code.NUCLEOTIDE_CHARACTERS_DOUBLE_RNA;
				}
				else if(character_grouping_number == 3) {
					return Code.NUCLEOTIDE_CHARACTERS_TRIPLE_RNA;
				}
				else {
					System.out.println("ERROR: Too high a grouping number! Exiting");
					System.exit(1);
					return null;
				}
			}
			else {
				System.out.println("ERROR: somehow ended up with an impossible case, exiting");
				System.exit(1);
				return null;
			}
		}
		else if(sequence_data_type == Code.AMINO_ACID) {
			return Code.AMINO_ACID_CHARACTERS_SINGLE;
		}
		else if(sequence_data_type == Code.STANDARD) {
			return Code.STANDARD_CHARACTERS_SINGLE;
		}
		else {
			System.out.println("ERROR: Un-handled sequence data type encountered in getMembers2(), exiting");
			System.exit(1);
			return null;
		}
	}
	
	
	/*
		Returns the total number of characters that make up a paricular code.
	*/
	public int getTotal()
	{
		if(sequence_data_type == Code.NUCLEOTIDE_DNA) {
			return (int)(Math.pow(Code.NUCLEOTIDE_CHARACTERS_SINGLE_DNA.length, character_grouping_number));
		}
		else if(sequence_data_type == Code.NUCLEOTIDE_RNA) {
			return (int)(Math.pow(Code.NUCLEOTIDE_CHARACTERS_SINGLE_RNA.length, character_grouping_number));
		}
		else if(sequence_data_type == Code.AMINO_ACID) {
			return (int)(Math.pow(Code.AMINO_ACID_CHARACTERS_SINGLE.length, character_grouping_number));
		}
		else {
			return members.size();
		}
	}
	
	
	/*
		uses instance variables for character grouping number and offset to determine the 
		grouped characters (members) present in a list of sequences, in the same manner
		the Alignment class creates sites lists of Sites. They could not be implemented
		simultaneously due to ordering issues
	*/
	public String[] determineContainedMembers(LinkedList<Sequence> sequences) {
		ListIterator<Sequence> iter = sequences.listIterator();
		int length = Alignment.getShortestSequenceLength(sequences); //length of shortest string in sequences
		int size = sequences.size(); //number of sequences in the list
		character_groups = new HashSet<String>();
		
		for(int i = offset_number; i < length; i += character_grouping_number) {
			if((i + character_grouping_number) < length)
			{
				iter = sequences.listIterator();
				for(int j = 0; j < size; j++)
				{
					String curr_string = ((Sequence)iter.next()).getSequence();
					String newBase = curr_string.substring(i, i+character_grouping_number);
					character_groups.add(newBase);
				}
			}
		}
		
		return character_groups.toArray(new String[] {});
	}
	
	
	/*
		Recodes a site at a particular position.
		Site s is made up of group characters.
		Character at position in the group will be recoded according to recode_index.
	*/
	public Site recode(Site s, int recode_index, int group, int position)
	{
		char[] replace_list;
		Site new_site = new Site(s.getPosition(), s.getSize());
		String temp = s.toString();
		
		// When sites are single nucleotides/amino acids
		if(group == 1)
		{
			// Handle all recoding schemes except those where two recoding steps are involved
			if(recode_index != 3 && recode_index != 6 && recode_index != 9)
			{
				replace_list = getReplacements(recode_index);
				
				for(int i = 1; i< replace_list.length; i++)
				{
					temp = temp.replace(replace_list[i],replace_list[0]);
				}
			} 
			// Handle recoding where two steps are involved.
			else 
			{
				replace_list = getReplacements(recode_index-2);
				for(int i = 1; i< replace_list.length; i++)
				{
					temp = temp.replace(replace_list[i],replace_list[0]);
				}
				replace_list = getReplacements(recode_index-1);
				for(int j = 1; j< replace_list.length; j++)
				{
					temp = temp.replace(replace_list[j],replace_list[0]);
				}
			}
			String character;
			for(int j = 0; j < s.getSize(); j++)
			{
				character = temp.substring(j,j+group);
				new_site.addBase(character);
			}
		}
		// When sites are di/tri nucleotides
		else
		{
			String pos = new String();
			for(int i = 0; i< temp.length();i+=group){
				pos+=temp.charAt(i+position);
			}
			if(recode_index != 3 && recode_index != 6 && recode_index != 9){
				replace_list = getReplacements(recode_index);
				for(int i = 1; i< replace_list.length; i++){
					pos = pos.replace(replace_list[i],replace_list[0]);
				}
			}
			else {
				replace_list = getReplacements(recode_index-2);
				for(int i = 1; i< replace_list.length; i++){
					pos = pos.replace(replace_list[i],replace_list[0]);
				}
				replace_list = getReplacements(recode_index-1);
				for(int j = 1; j< replace_list.length; j++){
					pos = pos.replace(replace_list[j],replace_list[0]);
				}
			}
			String character = new String();
			// When sites are dinucleotides
			if(group == 2){
				for(int k = 0; k < pos.length(); k++){
					if(position == 0){
						character = pos.substring(k, k+1) + temp.substring(2*k+1,2*k+2);
					} else {
						character = temp.substring(2*k,2*k+1) + pos.substring(k, k+1);
					}
					new_site.addBase(character);
				}
			}
			// When sites are trinucleotides
			else if(group == 3){
				for(int k = 0; k < pos.length(); k++){
					if(position == 0){
						character = pos.substring(k, k+1) + temp.substring(3*k+1,3*k+3);
					} else if(position == 1){
						character = temp.substring(3*k,3*k+1) + pos.substring(k, k+1) + temp.substring(3*k+2,3*k+3);
					}else{
						character = temp.substring(3*k,3*k+2) + pos.substring(k, k+1);
					}
					new_site.addBase(character);
				}
				System.out.println("");
			}
		}
		return new_site;
	}
	
	/*
		Returns the an array of the characters to be replaced and the character that will replace them
		The replacement character is the first in the array, while the remaining characters are the ones it is replacing
		Eg. For case 1, 'R' is replacing any 'G's or 'A's.
	*/
	public char[] getReplacements(int index)
	{
		char[] list = new char[4];
		switch(index)
		{
			case 0: break;
			case 1: list[0] = 'R'; list[1] ='G'; list[2] ='A'; break;
			case 2: list[0] = 'Y'; list[1] ='T'; list[2] ='C'; break;
			case 4: list[0] ='K'; list[1] ='G'; list[2]='T'; break;
			case 5: list[0] = 'M'; list[1]='A'; list[2]='C'; break;
			case 7: list[0] ='S'; list[1]='C'; list[2]='G'; break;
			case 8: list[0] = 'W'; list[1] ='A'; list[2] = 'T'; break;
			case 10: list[0] = 'B'; list[1] = 'C'; list[2] ='G'; list[3] = 'T'; break;
			case 11: list[0] = 'D'; list[1] = 'A'; list[2] ='G'; list[3] = 'T'; break;
			case 12: list[0] = 'H'; list[1] = 'A'; list[2] ='C'; list[3] = 'T'; break;
			case 13: list[0] = 'V'; list[1] = 'A'; list[2] ='C'; list[3] = 'G'; break;
		}
		return list;
	}
	
	
	/*
		Creates map of the ambiguous character ratio ie no of nucleotides each character represents.
		Used by Site:completeness method
	*/
	public Map<String,Integer> ambigRatio()
	{
		Map<String,Integer> m = new HashMap<String,Integer>();
		m.put("A",1);m.put("C",1);m.put("G",1);m.put("T",1);
		m.put("R",2);m.put("Y",2);m.put("K",2);m.put("M",2);m.put("S",2);m.put("W",2);
		m.put("B",3);m.put("D",3);m.put("H",3);m.put("V",3);
		return m;
	}
		
	
	/*
		Some helper 'get' methods
	*/
	public char getUnknown(){return unknown;}
	public char getMissing(){return missing;}
	public int getOffset() { return offset_number; }
	public int getGroupingNumber() { return character_grouping_number; }
}