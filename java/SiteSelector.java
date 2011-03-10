import java.awt.*;
import javax.swing.*;
import java.awt.event.*;

public class SiteSelector extends JPanel implements ActionListener {
	
	private int[] cp_numbers;
	private int[] mnic_numbers;
	private double[] entropy_numbers;
	private int[] compl_numbers;
	private JPanel radioPanel;
	private JRadioButton cpOne, cpTwo, cpThree, complete, incomplete;
	//private JTextField mnic_min, mnic_max, ent_min, ent_max;
	//private JButton set_button;
	private RecodeSelector recode;
	private boolean altered;
	
	public SiteSelector(RecodeSelector rs) {
		super(new BorderLayout());
		cp_numbers = new int[3];
		//mnic_numbers = new int[2];
		//entropy_numbers = new double[2];
		compl_numbers = new int[2];
		cp_numbers[0]=1;cp_numbers[1]=1;cp_numbers[2]=1;
		//mnic_numbers[0]=1;mnic_numbers[1]=19;
		//entropy_numbers[0]=0;mnic_numbers[1]=100;
		compl_numbers[0]=1;compl_numbers[1]=1;
		recode = rs;
		altered = false;
		
		/*int[] mnic_options = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
		JComboBox mnic_list = new JComboBox(mnic_options);
		mnic_list.setEditable(true);
		mnic_list.addActionListener(this);*/
		
		//creates radio buttons for codon position
				
		cpOne = new JRadioButton("1st");
		cpOne.setActionCommand("0");
		cpOne.setSelected(true);
		
		cpTwo = new JRadioButton("2nd");
		cpTwo.setActionCommand("1");
		cpTwo.setSelected(true);
		
		cpThree = new JRadioButton("3rd");
		cpThree.setActionCommand("2");
		cpThree.setSelected(true);

		//mnic_min = new JTextField("1");
		//mnic_min.setEnabled(true);
		//mnic_max = new JTextField("19");
		//mnic_max.setEnabled(true);
		//ent_min = new JTextField("0");
		//ent_min.setEnabled(true);
		//ent_max = new JTextField("100");
		//ent_max.setEnabled(true);
		//set_button=new JButton("Set");
		
		complete = new JRadioButton("Complete");
		complete.setActionCommand("0");
		complete.setSelected(true);
		
		incomplete = new JRadioButton("Incomplete");
		incomplete.setActionCommand("1");
		incomplete.setSelected(true);
		
				
		//set ActionListeners for all of the bin grouping radiobuttons
		cpOne.addActionListener(this);
		cpTwo.addActionListener(this);
		cpThree.addActionListener(this);
		//set_button.addActionListener(new SetListener());
		complete.addActionListener(new ComplListener());
		incomplete.addActionListener(new ComplListener());
						
		//puts the radio buttons in a column in a panel
		radioPanel = new JPanel(new GridLayout(12,2));
		radioPanel.add(new JLabel("Codon Position:"));
		radioPanel.add(new JLabel(""));
		radioPanel.add(cpOne);
		radioPanel.add(new JLabel(""));
		radioPanel.add(cpTwo);
		radioPanel.add(new JLabel(""));
		radioPanel.add(cpThree);
		radioPanel.add(new JLabel(""));
		//radioPanel.add(new JLabel("MNIC:"));
		//radioPanel.add(new JLabel(""));
		//radioPanel.add(new JLabel("Min"));
		//radioPanel.add(new JLabel("Max"));
		//radioPanel.add(mnic_min);
		//radioPanel.add(mnic_max);
		//radioPanel.add(new JLabel("Entropy:"));
		//radioPanel.add(new JLabel(""));
		//radioPanel.add(new JLabel("Min"));
		//radioPanel.add(new JLabel("Max"));
		//radioPanel.add(ent_min);
		//radioPanel.add(ent_max);
		radioPanel.add(new JLabel("Completeness:"));
		radioPanel.add(new JLabel(""));
		radioPanel.add(complete);
		radioPanel.add(new JLabel(""));
		radioPanel.add(incomplete);
		radioPanel.add(new JLabel(""));
		radioPanel.add(new JLabel(""));
		radioPanel.add(new JLabel(""));
		radioPanel.add(new JLabel(""));
		radioPanel.add(new JLabel(""));
		radioPanel.add(new JLabel(""));
		radioPanel.add(new JLabel(""));
		radioPanel.add(new JLabel(""));
		radioPanel.add(new JLabel(""));
		radioPanel.add(new JLabel(""));
		radioPanel.add(new JLabel(""));

		//radioPanel.add(new JLabel(""));
		//radioPanel.add(set_button);
		
		//adds radioPanel to the BinSelector
		add(radioPanel, BorderLayout.CENTER);
	}
	
	public void actionPerformed(ActionEvent e)
	{
		altered=true;
		if(cp_numbers[Integer.parseInt(e.getActionCommand())]==1){
			cp_numbers[Integer.parseInt(e.getActionCommand())]=0;
		}else{
			cp_numbers[Integer.parseInt(e.getActionCommand())]=1;
		}
	}
	
	//create ActionListeners for set button
	/*class SetListener implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			altered=true;
			String str = null;
			try{str = mnic_min.getText();}
			catch(NullPointerException npe){}
			if(str!=null){mnic_numbers[0] = new Integer(str.trim());}
			try{str = mnic_max.getText();}
			catch(NullPointerException npe){}
			if(str!=null){mnic_numbers[1] = new Integer(str.trim());}
			//try{str = ent_min.getText();}
			//catch(NullPointerException npe){}
			//if(str!=null){entropy_numbers[0] = new Double(str.trim());}
			//try{str = ent_max.getText();}
			//catch(NullPointerException npe){}
			//if(str!=null){entropy_numbers[1] = new Double(str.trim());}
		}
	}*/
	
	//create ActionListeners for completeness
	class ComplListener implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			altered=true;
			if(compl_numbers[Integer.parseInt(e.getActionCommand())]==1){
				compl_numbers[Integer.parseInt(e.getActionCommand())]=0;
			}else{
				compl_numbers[Integer.parseInt(e.getActionCommand())]=1;
			}
		}
	}
	
	public void updatePos(int pos)
	{
		if(pos==0){
			if(cpOne.isSelected()){
				cp_numbers[0]=0;
			}else{
				cp_numbers[0]=1;
			}
		}
		else if(pos==1){
			if(cpTwo.isSelected()){
				cp_numbers[1]=0;
			}else{
				cp_numbers[1]=1;
			}
		}
		else if(pos==2){
			if(cpThree.isSelected()){
				cp_numbers[2]=0;
			}else{
				cp_numbers[2]=1;
			}
		}
	}
	
	public void updateState(int code_id, int grouping)
	{
		System.out.println("Grouping = " + grouping);
		if(code_id == 1 && grouping == 1){
			cpOne.setEnabled(true);
			cpTwo.setEnabled(true);
			cpThree.setEnabled(true);
		}
		else if(code_id == 2 && grouping == 1){
			cpOne.setEnabled(true);
			cpTwo.setEnabled(true);
			cpThree.setEnabled(false);
		}else{
			cpOne.setSelected(true);
			cpTwo.setSelected(true);
			cpThree.setSelected(true);
			cpOne.setEnabled(false);
			cpTwo.setEnabled(false);
			cpThree.setEnabled(false);
		}
	}
	
	public int[] getcpNumbers() {return cp_numbers; }	
	//public int[] getmnicNumbers() { return mnic_numbers; }
	//public double[] getentropyNumbers() { return entropy_numbers; }
	public int[] getcomplNumbers() { return compl_numbers; }
	public boolean isAltered(){ return altered; }
}