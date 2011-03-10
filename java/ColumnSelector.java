import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

/*
	Column Selector is the GUI panel where users decide how they want to order the sites
	in the compatibility plot
*/
public class ColumnSelector extends JPanel implements ActionListener
{
	private int current_selection;
	public static final int AS_IS = 0;
	public static final int REORDER_CHANGES = 1;
	public static final int REORDER_CHANGES_SPLITS = 2;
	public static final int REORDER_COMP = 3;
	
	/* 
		Class Declaration
	*/
	public ColumnSelector()
	{
		super(new BorderLayout());
		current_selection = ColumnSelector.AS_IS;
		
		//create radiobuttons
		JRadioButton as_is_button = new JRadioButton("Leave as is");
		as_is_button.setActionCommand(Integer.toString(ColumnSelector.AS_IS));
		as_is_button.setSelected(true);
		
		JRadioButton reorder_changes_button = new JRadioButton("MNIC (1 -> n)");
		reorder_changes_button.setActionCommand(Integer.toString(ColumnSelector.REORDER_CHANGES));
		reorder_changes_button.setSelected(true);
		
		JRadioButton reorder_splits_button = new JRadioButton("Entropy (0 -> m)");
		reorder_splits_button.setActionCommand(Integer.toString(ColumnSelector.REORDER_CHANGES_SPLITS));
		reorder_splits_button.setSelected(true);
		
		JRadioButton reorder_comp_button = new JRadioButton("Compatibility (1 -> 0)");
		reorder_comp_button.setActionCommand(Integer.toString(ColumnSelector.REORDER_COMP));
		reorder_comp_button.setSelected(true);
		
		//puts the buttons into a ButtonGroup
		ButtonGroup col_group = new ButtonGroup();
		col_group.add(as_is_button);
		col_group.add(reorder_changes_button);
		col_group.add(reorder_splits_button);
		col_group.add(reorder_comp_button);
		
		//adds actionListeners for the buttons
		as_is_button.addActionListener(this);
		reorder_changes_button.addActionListener(this);
		reorder_splits_button.addActionListener(this);
		reorder_comp_button.addActionListener(this);
		
		//adds buttons to a panel in a column
		JPanel radioPanel = new JPanel(new GridLayout(12,1));
		radioPanel.add(as_is_button);
		radioPanel.add(reorder_changes_button);
		radioPanel.add(reorder_splits_button);
		radioPanel.add(reorder_comp_button);
		
		//adds radioPanel to this RecodeSelector
		add(radioPanel, BorderLayout.CENTER);
		
	}
	
	/*
		Main action listener
		Is called when the user changes the radio button in this panel of the GUI
	*/
	public void actionPerformed(ActionEvent e) 
	{
		int selection = Integer.parseInt(e.getActionCommand());
		current_selection = selection;
	}
	
	/*
		Returns the number of the current selection of how the sites are to be ordered.
	*/
	public int getColumnSelection() {
		return current_selection;
	}
}