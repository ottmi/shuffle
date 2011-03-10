import java.awt.*;
import java.awt.event.*;
import java.util.*;
import javax.swing.*;
import java.lang.*;
import java.io.*;

/*
	This class brings together and integrates all the separate components
	of the GUI, laying them out in a JFrame which is intended for use as 
	the main Frame of a program
*/
public class ApplicationGui extends JFrame {
	
	GridTabbedPane tp_grids;
	JTabbedPane tp_selectors;
	SequenceSelector2 ss;
	StatsPane2 stats_pane;
	BinSelector binSelector;
	RecodeSelector recodeSelector;
	ColumnSelector colSelector;
	OutputSelector outSelector;
	SiteSelector siteSelector;
	private JFileChooser fileChooser;
	int grid_counter;
	private Alignment al;
	private String format;
	private File curr_file;
	private Calendar cal;
	private int pixelSize;
	
	/* 
		Contructor: Creates a new instance based on the name of the frame - "Shuffle"
	*/
	public ApplicationGui(String frameName) 
	{
		super(frameName);
		Container contentPane = getContentPane();
		
		
		makeFileChooser(); //sets up JFileChooser for opening files
		
		JFrame appFrame = new JFrame("Shuffle"); //the main frame for the GUI
		
		// Initialise all the selectors and dependant variables
		recodeSelector = new RecodeSelector();
		siteSelector = new SiteSelector(recodeSelector);
		binSelector = new BinSelector(recodeSelector,siteSelector);
		colSelector = new ColumnSelector();
		outSelector = new OutputSelector();
		ss = new SequenceSelector2();
		stats_pane = new StatsPane2();
		cal = Calendar.getInstance();
		grid_counter = 0;
		pixelSize = 1;
		
		/*
			Listener for the draw plot button/menuItem
		*/
		class GOListener implements ActionListener 
		{
			public void actionPerformed(ActionEvent event) 
			{
				BinaryGridImage new_grid = GO();
				grid_counter++;
				stats_pane.appendln("Run " + grid_counter);
				if(new_grid != null) {
					tp_grids.addTab("Grid " + grid_counter , new JScrollPane(new_grid));
					stats_pane.appendln("Plot generated");
				}
				else{ stats_pane.appendln("No plot generated");}
				
				stats_pane.appendln("	File: " + fileChooser.getName(curr_file));
				stats_pane.appendln("	Selected Informative sites: " + al.selectedTotal());
				stats_pane.appendln("	Overall Compatibility: " + new_grid.getOverallComp());		
									
				if(binSelector.getBinSelectionNumber() > 0 && binSelector.getBinSelectionNumber() < 4){
					stats_pane.append("	Recode: ");
					if(binSelector.getBinSelectionNumber() == 1){
						stats_pane.appendln("(" + recodeSelector.getSelectionName(0) + ")");
					}
					if(binSelector.getBinSelectionNumber() == 2){
						stats_pane.append("Position 1: (" + recodeSelector.getSelectionName(0));
						stats_pane.appendln(") Position 2: (" + recodeSelector.getSelectionName(1) + ")");
					}
					if(binSelector.getBinSelectionNumber() == 3){
						stats_pane.append("Position 1: (" + recodeSelector.getSelectionName(0));
						stats_pane.append(") Position 2: (" + recodeSelector.getSelectionName(1));
						stats_pane.appendln(") Position 3: (" + recodeSelector.getSelectionName(2) + ")");
					}
				}
				if(new_grid != null) {	
					if(colSelector.getColumnSelection() == 0){
						stats_pane.appendln("	Site's order: As is");
					}else if(colSelector.getColumnSelection() == 1){
						stats_pane.appendln("	Site's order: MNIC(0->n)");
					}else if(colSelector.getColumnSelection() == 2){
						stats_pane.appendln("	Site's order: Compatibility(1->0)");
					}
				}
					
				if(!siteSelector.isAltered()){
					stats_pane.appendln("	All sites selected");
				}else{
					stats_pane.append("	Selected Sites: ");
					int[] cps = siteSelector.getcpNumbers();
					for(int i=0; i<3; i++){
						if(cps[i] == 1){ stats_pane.append ("CP" + (i+1) + " ");}
					}
					int[] compl = siteSelector.getcomplNumbers();
					if(compl[0] == 1){ stats_pane.append("complete ");}	
					if(compl[1] == 1){ stats_pane.append("incomplete ");}
					stats_pane.appendln("");
				}
				new_grid = null;
			}
		}
		
		tp_grids = new GridTabbedPane(new GOListener()); //this will be loaded into the app frame
		tp_grids.setPreferredSize(new Dimension(400,400));
		
		
		//sets up tp_selectors, which is a tabbed pane containing all the selector GUI components
		tp_selectors = new JTabbedPane();
		tp_selectors.addTab("Taxa", ss);
		tp_selectors.addTab("Data Type", binSelector);
		tp_selectors.addTab("Recode", recodeSelector);
		tp_selectors.addTab("Site Selection", siteSelector);
		tp_selectors.addTab("Order of Sites", colSelector);
		tp_selectors.addTab("Output File", outSelector);
		
		//contentPane loading using GridBagLayout
		contentPane.setLayout(new GridBagLayout());
		GridBagConstraints constraints;
		
		//add sequence selector at top left
		//sets position to top left
		constraints = new GridBagConstraints();
		constraints.gridx = 0;
		constraints.gridy = 0;
		constraints.weighty = 0;	//sets weight for vertical stretching
		constraints.fill = GridBagConstraints.VERTICAL;
		contentPane.add(tp_selectors, constraints);
		
		//add gridtabbedpane to top right
		constraints = new GridBagConstraints();
		constraints.gridx = 1;
		constraints.gridy = 0;
		constraints.weightx = 0.5;
		constraints.weighty = 0;
		constraints.fill = GridBagConstraints.BOTH;
		contentPane.add(tp_grids, constraints);
		
		//add statspane along bottom
		constraints = new GridBagConstraints();
		constraints.gridx = 0;
		constraints.gridy = 1;
		constraints.gridwidth = 2;	//makes component stretch across two cells
		constraints.ipady = 50;
		constraints.weightx = 0.5;
		constraints.weighty = 0;
		constraints.fill = GridBagConstraints.BOTH;
		contentPane.add(stats_pane, constraints);
		
		
		
		/*
			Listener of the View tab in the menu bar
		*/
		class ViewMenuListener implements ActionListener {
			public void actionPerformed(ActionEvent event) {
				String command = event.getActionCommand();
//				System.out.println(command);
				selectSelector(command);
			}
		}
		
		//Assemble the JMenuBar
		JMenuBar menuBar = new JMenuBar();
		
		JMenu fileMenu = new JMenu("File");
		fileMenu.setMnemonic(KeyEvent.VK_F);
		menuBar.add(fileMenu);
		
		JMenuItem openFile = new JMenuItem("Open File", KeyEvent.VK_O);
		openFile.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_O, ActionEvent.ALT_MASK));
		openFile.addActionListener(new OpenFileListener());
		fileMenu.add(openFile);
		
		JMenuItem exitProgram = new JMenuItem("Exit", KeyEvent.VK_Q);
		exitProgram.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_Q, ActionEvent.META_MASK));
		exitProgram.addActionListener(new ExitProgramListener());
		fileMenu.add(exitProgram);
		
		JMenu viewMenu = new JMenu("View");
		viewMenu.setMnemonic(KeyEvent.VK_V);
		menuBar.add(viewMenu);
		
		JMenuItem taxaItem = new JMenuItem("Taxa", KeyEvent.VK_T);
		taxaItem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_T, ActionEvent.META_MASK));
		taxaItem.setActionCommand("Taxa");
		taxaItem.addActionListener(new ViewMenuListener());
		viewMenu.add(taxaItem);
		
		JMenuItem dataItem = new JMenuItem("Data", KeyEvent.VK_D);
		dataItem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_D, ActionEvent.META_MASK));
		dataItem.setActionCommand("Data Type");
		dataItem.addActionListener(new ViewMenuListener());
		viewMenu.add(dataItem);
		
		JMenuItem siteItem = new JMenuItem("Site Selector", KeyEvent.VK_C);
		siteItem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_C, ActionEvent.META_MASK));
		siteItem.setActionCommand("Site Selector");
		siteItem.addActionListener(new ViewMenuListener());
		viewMenu.add(siteItem);
		
		JMenuItem recodeItem = new JMenuItem("Recode", KeyEvent.VK_R);
		recodeItem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_R, ActionEvent.META_MASK));
		recodeItem.setActionCommand("Recode");
		recodeItem.addActionListener(new ViewMenuListener());
		viewMenu.add(recodeItem);
		
		JMenuItem orderItem = new JMenuItem("Order of Sites", KeyEvent.VK_S);
		orderItem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_S, ActionEvent.META_MASK));
		orderItem.setActionCommand("Order of Sites");
		orderItem.addActionListener(new ViewMenuListener());
		viewMenu.add(orderItem);
		
		JMenuItem outputItem = new JMenuItem("Output File", KeyEvent.VK_O);
		outputItem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_O, ActionEvent.META_MASK));
		outputItem.setActionCommand("Output File");
		outputItem.addActionListener(new ViewMenuListener());
		viewMenu.add(outputItem);
		
		JMenu gridMenu = new JMenu("Plot");
		gridMenu.setMnemonic(KeyEvent.VK_G);
		menuBar.add(gridMenu);
		
		JMenuItem drawNewGrid = new JMenuItem("Draw Plot", KeyEvent.VK_N);
		drawNewGrid.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_N, ActionEvent.ALT_MASK));
		drawNewGrid.addActionListener(new GOListener());
		gridMenu.add(drawNewGrid);
		
		JMenu pixelSizeMenu = new JMenu("Pixel Size");
		JRadioButtonMenuItem pixel1 = new JRadioButtonMenuItem("1");
		JRadioButtonMenuItem pixel2 = new JRadioButtonMenuItem("2");
		JRadioButtonMenuItem pixel3 = new JRadioButtonMenuItem("3");
		JRadioButtonMenuItem pixel4 = new JRadioButtonMenuItem("4");
		JRadioButtonMenuItem pixel5 = new JRadioButtonMenuItem("5");
		JRadioButtonMenuItem pixel10 = new JRadioButtonMenuItem("10");
		JRadioButtonMenuItem pixel15 = new JRadioButtonMenuItem("15");
		JRadioButtonMenuItem pixel20 = new JRadioButtonMenuItem("20");
		JRadioButtonMenuItem pixel30 = new JRadioButtonMenuItem("30");
		
		pixel1.setSelected(true);
		
		pixel1.setActionCommand("1");
		pixel2.setActionCommand("2");
		pixel3.setActionCommand("3");
		pixel4.setActionCommand("4");
		pixel5.setActionCommand("5");
		pixel10.setActionCommand("10");
		pixel15.setActionCommand("15");
		pixel20.setActionCommand("20");
		pixel30.setActionCommand("30");
		
		pixelSizeMenu.add(pixel1);
		pixelSizeMenu.add(pixel2);
		pixelSizeMenu.add(pixel3);
		pixelSizeMenu.add(pixel4);
		pixelSizeMenu.add(pixel5);
		pixelSizeMenu.add(pixel10);
		pixelSizeMenu.add(pixel15);
		pixelSizeMenu.add(pixel20);
		pixelSizeMenu.add(pixel30);
		
		gridMenu.add(pixelSizeMenu);
		
		ButtonGroup pixelGroup = new ButtonGroup();
		pixelGroup.add(pixel1);
		pixelGroup.add(pixel2);
		pixelGroup.add(pixel3);
		pixelGroup.add(pixel4);
		pixelGroup.add(pixel5);
		pixelGroup.add(pixel10);
		pixelGroup.add(pixel15);
		pixelGroup.add(pixel20);
		pixelGroup.add(pixel30);
		
		PixelSizeListener psl = new PixelSizeListener();
		pixel1.addActionListener(psl);
		pixel2.addActionListener(psl);
		pixel3.addActionListener(psl);
		pixel4.addActionListener(psl);
		pixel5.addActionListener(psl);
		pixel10.addActionListener(psl);
		pixel15.addActionListener(psl);
		pixel20.addActionListener(psl);
		pixel30.addActionListener(psl);

		setJMenuBar(menuBar);
	}
	
	/*
		Sets up the open file JFileChooser as desired
	*/
	private void makeFileChooser() {
		fileChooser = new JFileChooser();
		fileChooser.setAcceptAllFileFilterUsed(false);
		fileChooser.addChoosableFileFilter(new NexusFilter());
		fileChooser.addChoosableFileFilter(new InterleavedPhylipFilter());
		fileChooser.addChoosableFileFilter(new SequentialPhylipFilter());
		fileChooser.addChoosableFileFilter(new FastAFilter());
	}

/*
		Helper method for selecting the appropriate Selector component within tp_selectors
	*/
	private void selectSelector(String title) {
		tp_selectors.setSelectedIndex(tp_selectors.indexOfTab(title));
	}
	
	/*
		Listens for when the user wishes to open a new file
	*/
	class OpenFileListener implements ActionListener {
		public void actionPerformed(ActionEvent e) {
//			System.out.println("OpenFileListener triggered");
			int returnVal = fileChooser.showOpenDialog(ApplicationGui.this);
			if(returnVal == JFileChooser.APPROVE_OPTION) {
				File file = fileChooser.getSelectedFile();
				curr_file = file;
				javax.swing.filechooser.FileFilter filter = fileChooser.getFileFilter();
				createNewSequenceSelector(file, filter.getClass());
			}
		}
	}
	
	
	/*
		Listens for when the user wishes to close the program
		This was necessary as sometimes OS specific close key combinations (such as
		command + q in Mac OS X) would not be recognised
	*/
	class ExitProgramListener implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			System.exit(0);
		}
	}
	
	
	/*
		Listens for when the user selects a different number of pixels to
		be the width of boxes in the plots drawn
	*/
	class PixelSizeListener implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			int selection = Integer.parseInt(e.getActionCommand());
			setPixelSize(selection);
		}
	}
	//helper method
	private void setPixelSize(int size) {
		pixelSize = size;
	}

	/*
		Creates a new SequenceSelector object of the appropriate file format type.
		Uses JFileChooser filter to determine which parser to use, then uses
		the parsed sequences to construct a new SequenceSelector
	*/
	private void createNewSequenceSelector(File file, Class class_type) 
	{
		LinkedList<Sequence> seqs = new LinkedList<Sequence>();
		if(class_type == NexusFilter.class) {
			NexusParser nex = new NexusParser(file);
			try{
				seqs = nex.parseFile();
			} 
			catch(Exception ex)
			{
				stats_pane.appendln("File is not in the nexus format. Please try again");
			}
			format = "NEXUS";
		}
		else if(class_type == SequentialPhylipFilter.class) {
			
			PhylipParser phy = new PhylipParser(file,'S');
			try{
				seqs = phy.parseFile();
			} 
			catch(Exception ex)
			{
				stats_pane.appendln("File is not in the sequential phylip format. Please try again");
			}
			format = "SPHYLIP";
		}
		else if(class_type == InterleavedPhylipFilter.class) {
			PhylipParser phy = new PhylipParser(file,'I');
			try{
				seqs = phy.parseFile();
			} 
			catch(Exception ex)
			{
				stats_pane.appendln("File is not in the interleaved phylip format. Please try again");
			}
			format = "IPHYLIP";
		}
		else if(class_type == FastAFilter.class) {
			FastAParser fap = new FastAParser(file);
			try{
				seqs = fap.parseFile();
			} 
			catch(Exception ex)
			{
				stats_pane.appendln("File is not in the FASTA format. Please try again");
			}
			format = "FASTA";
		}
		else {
			System.exit(1);
		}

		tp_selectors.remove(ss);
		ss = new SequenceSelector2(seqs);
		tp_selectors.add(ss, 0);
		tp_selectors.setTitleAt(0, "Taxa");
		tp_selectors.setSelectedIndex(0);
	}
					
					

	/*
		Method which brings together information from all the selectors and makes the necessary calculations
	*/
	private BinaryGridImage GO() 
	{
		LinkedList<Sequence> select = ss.getSelectedSequences();
		if(select.size() == 0) {
			stats_pane.appendln("no sequences selected!");
			return null;
		}
//		System.out.println("Selected Sequences = " + select.size());
		int dataType = binSelector.getBinSelectionNumber();
//		System.out.println("Data type = " + dataType);
		int codeType;
		int maxMNIC = 0;
		char uk; 
		char mis = '-';
		int grouping = 1;
		int frame = 0;
		
		//determine type of sequence data
		if(dataType == 0){
			codeType = Code.AMINO_ACID;
			uk = 'X';
			maxMNIC = 21;
		}
		else if(dataType > 0 && dataType < 4){
			codeType = Code.NUCLEOTIDE_DNA;
			uk = 'N';
			if(binSelector.getGrouping() == 0){
				grouping  = binSelector.getBinSelectionNumber();
			}
			frame = binSelector.getBinOffsetNumber();
			maxMNIC = 3;
		}
		else if(dataType==4){
			codeType = Code.STANDARD;
			uk='?';
			maxMNIC=26+10-1;
		}
		else{
			codeType = Code.NUCLEOTIDE_RNA;
			uk = 'N';
			if(binSelector.getGrouping() == 0){
				grouping  = binSelector.getBinSelectionNumber();
			}
			frame = binSelector.getBinOffsetNumber();
			maxMNIC = 3;
		}
		
		System.out.println("Grouping = " + grouping);
		
		Code new_code = new Code("A", codeType,uk, mis, grouping, frame); //Generate new Code object reflecting the sequence type
		//String[] determinedMembers = new_code.determineContainedMembers(select); //find all possible site configurations
		al = new Alignment(select,new_code,recodeSelector); //turn selected sequences into an alignment
		LinkedList<Site> sites = al.findInformative(); //find the informative sites in the alignment
		
		System.out.println("Informative Total = " + al.informativeTotal());
	
		// If only certain sites are to be analysed...
		LinkedList<Site> sel_sites = new LinkedList<Site>();
		if(siteSelector.isAltered()){
			int[] cp = siteSelector.getcpNumbers();
			int[] compl = siteSelector.getcomplNumbers();
			ListIterator<Site> iter = sites.listIterator();
			while(iter.hasNext()){
				Site site_cur = iter.next();
				int site_cp = site_cur.getCP(binSelector.getBinSelectionNumber());
				int site_mnic = site_cur.getChanges();
				boolean site_compl = site_cur.isComplete();
				if((site_compl && compl[0]==1) || (!site_compl && compl[1]==1)){
					if(grouping == 1){
						for(int i=0;i<cp.length;i++){
							if((cp[i]==1 && site_cp == i+1)){
								sel_sites.add(site_cur);
								break;
							}
						}
					}else{sel_sites.add(site_cur);}
				}
			}
		}else{sel_sites = sites;}
		al.setSelected(sel_sites);
		
		//double nss = null;
		
		BinaryGridImage bgi = null;
		
		// generates the boolean 2D array used to generate compatibility plot
		if(outSelector.getDraw()){
			System.out.print("Drawing Plot...");
			//orders sites by user selected criterion
			LinkedList<Site> osites = new LinkedList<Site>();
			if(colSelector.getColumnSelection() == 1){
				osites = al.orderSitesByChanges();
			}
			else if(colSelector.getColumnSelection() == 2){
				osites = al.orderSitesByEntropy();
			}
			if(colSelector.getColumnSelection() == 3){
				osites = al.orderSitesByComp();
			}
			else if(colSelector.getColumnSelection() == 0){
				osites = sel_sites;
			}
			boolean[][] plot = al.plot(osites, osites);
			try{
				bgi = new BinaryGridImage(plot, pixelSize, osites, osites, stats_pane);
			}
			catch(java.lang.OutOfMemoryError e) {
				stats_pane.appendln("Ran out of memory trying to draw grid");
				e.printStackTrace();
			}
			System.out.println("Done.");
		}
		else{
			al.calcComp();
		}
		
		//handles file output if selected
		if(outSelector.getOutput()){
			System.out.print("Outputting enhanced alignment...");
			LinkedList<Double> POCthresholds = outSelector.getPOCThreshs();
			double mincompthresh = outSelector.getCompThresh();
			int changes = outSelector.getMaxChanges();
			if(changes > maxMNIC){
				changes = maxMNIC;
			}
			int gens = outSelector.getGenerations();
			al.outputFile("SPHYLIP",mincompthresh,POCthresholds,gens,changes,outSelector.getFilename());
			outSelector.resetSiglevels();
			System.out.println("Done.");
		}
		
		/*if(outSelector.getNSS()){
			nss = al.NeighbourSimScore(plot, 1000);
			stats_pane.appendln("	NeighbourSimScore: " + nss);
		}*/
		
		if(outSelector.getSummary()){
			System.out.print("Printing site summary file...");
			al.printSiteSummary();
			System.out.println("Done.");
		}
		return bgi;
	}
}
			
			
			
			