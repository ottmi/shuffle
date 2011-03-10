import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.imageio.*;
import java.lang.Object;
import java.io.*;
import javax.swing.filechooser.*;

/*
	Class which is a component containing a JTabbedPane for showing BinaryGridImages, and a JPanel with buttons for activities related to the images
*/
public class GridTabbedPane extends JComponent {
	
	private JTabbedPane tabs;
	private JPanel buttonPanel;
	private int grid_counter;
	private final JFileChooser fileChooser;
	
	public GridTabbedPane(ActionListener drawListener) {
		tabs = new JTabbedPane();
		grid_counter = 0;
		
		//instantiates and sets up JFileChooser for saving images
		fileChooser = new JFileChooser();
		fileChooser.setAcceptAllFileFilterUsed(false);
		JPGFilter jpg = new JPGFilter();
		fileChooser.addChoosableFileFilter(jpg);
		fileChooser.addChoosableFileFilter(new PNGFilter());
		fileChooser.setFileFilter(jpg);
		fileChooser.setApproveButtonText("Save");
		
		//create ActionListeners for button panel
		class CloseTabListener implements ActionListener {
			public void actionPerformed(ActionEvent e) {
				System.out.println("close tablistener invoked");
				if(tabs.getTabCount() > 0) {
					tabs.remove(tabs.getSelectedIndex());
				}
			}
		}
		
		/*
			Performs actions when 'Save Image' button activated
		*/
		class SaveImageListener implements ActionListener 
		{
			public void actionPerformed(ActionEvent e) 
			{
				if(tabs.getTabCount() > 0) {
					
					//retrieve the BinaryGridImage in the currently selected tab
					JViewport curr_component = (JViewport)((JScrollPane)tabs.getSelectedComponent()).getViewport();
					BinaryGridImage original_image = (BinaryGridImage)curr_component.getView();
					
					//display JFileChooser dialog
					int returnVal = fileChooser.showSaveDialog(GridTabbedPane.this);
					if(returnVal == JFileChooser.APPROVE_OPTION) {
						File f = fileChooser.getSelectedFile(); //filename and directory location where the user wishes the file go
						javax.swing.filechooser.FileFilter filter = fileChooser.getFileFilter(); //the file type they selected
						
						//checks to see if user has added the appropriate file extension already
						//first get the string after the last '.', as it could be a file extension
						String ext = null;
						String filename = f.getName();
//						System.out.println(filename);
						int i = filename.lastIndexOf('.');
						if(i > 0 && i < filename.length() - 1) {
							ext = filename.substring(i+1).toLowerCase();
						}
						
						//next, see if its the right extension for the selected filter, if not, save!
						if(filter.getClass() == PNGFilter.class) {
//							System.out.println("its a png!");
							if(ext != null && ext.equals(PNGFilter.getExtensionString())) { //user supplied correct extension
								original_image.writeImageToFile(f, PNGFilter.getExtensionString());
							}
							else {
								original_image.writeImageToFile(new File(f.getAbsolutePath()+"."+PNGFilter.getExtensionString()), PNGFilter.getExtensionString());
							}
						}
						else if(filter.getClass() == JPGFilter.class) {
//							System.out.println("its a jpg!");
							if(ext != null && ext.equals(JPGFilter.getExtensionString())) { //user supplied correct extension
								original_image.writeImageToFile(f, JPGFilter.getExtensionString());
							}
							else
								original_image.writeImageToFile(new File(f.getAbsolutePath()+"."+JPGFilter.getExtensionString()), JPGFilter.getExtensionString());
						}
						else if(filter.getClass() == GIFFilter.class) {
//							System.out.println("its a gif!");
							if(ext != null && ext.equals(GIFFilter.getExtensionString())) { //user supplied correct extension
								original_image.writeImageToFile(f, GIFFilter.getExtensionString());
							}
							else
								original_image.writeImageToFile(new File(f.getAbsolutePath()+"."+GIFFilter.getExtensionString()), GIFFilter.getExtensionString());
						}
						else { //unreachable!
							System.out.println("oh noes! this shouldn't be printing!");
						}
					}
				}
			}
		}
		
		//create a small class for the buttonpanel
		class GridTabbedPaneButtonPanel extends JPanel {
			public GridTabbedPaneButtonPanel(ActionListener draw_listener, ActionListener close_listener, ActionListener save_listener) {
				JButton draw = new JButton("Go!");
				JButton close = new JButton("Close Tab");
				JButton saveImage = new JButton("Save Image");
				
				//add ActionListeners to the buttons
				draw.addActionListener(draw_listener);
				close.addActionListener(close_listener);
				saveImage.addActionListener(save_listener);
				
				//layout and load
				setLayout(new GridLayout(1,0));
				add(draw);
				add(saveImage);
				add(close);
			}
		}
		
		//now that those thingies have been declared, we can instantiate a button panel
		buttonPanel = new GridTabbedPaneButtonPanel(drawListener, new CloseTabListener(), new SaveImageListener());
		
		//layout and load this GridTabbedPane object
		setLayout(new BorderLayout());
		add(tabs, BorderLayout.CENTER);
		add(buttonPanel, BorderLayout.SOUTH);
	}
	
	
	/*
		Adds a tab to the JTabbedPane, selects to the added tab
	*/
	public void addTab(String title, Component component) {
		tabs.addTab(title, component);
		tabs.setSelectedComponent(component);
	}
}

/*
	Filters for the JFileChooser dialog
*/
class JPGFilter extends javax.swing.filechooser.FileFilter {
	public boolean accept(File f) {
		return true;
	}
	public String getDescription() {
		return "JPEG (*.jpg)";
	}
	public static String getExtensionString() {
		return "jpg";
	}
}
class PNGFilter extends javax.swing.filechooser.FileFilter {
	public boolean accept(File f) {
		return true;
	}
	public String getDescription() {
		return "PNG (*.png)";
	}
	public static String getExtensionString() {
		return "png";
	}
}
class GIFFilter extends javax.swing.filechooser.FileFilter {
	public boolean accept(File f) {
		return true;
	}
	public String getDescription() {
		return "GIF (*.gif)";
	}
	public static String getExtensionString() {
		return "gif";
	}
}