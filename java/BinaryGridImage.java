import java.awt.*;
import javax.swing.*;
import java.awt.geom.Rectangle2D;
import java.awt.event.*;
import java.awt.image.*;
import java.io.*;
import javax.imageio.*;
import java.util.LinkedList;

/*
	A panel that displays a representation of a boolean array as white and black squares
*/
public class BinaryGridImage extends JPanel {
	
	static final long serialVersionUID = 1L;
	boolean[][] data;
	int width, height;
	int xmax, xmin, ymax,ymin;
	BufferedImage bimage;
	private final Site[] x_data, y_data;
	double overall_comp;
	
	/*
		constructs a new BinaryGridImage the size of which is the dimensions of the supplied array multiplied
		the square width, in pixels
	*/
	public BinaryGridImage(boolean[][] data, final int squareWidth, LinkedList<Site> x_datalist, LinkedList<Site> y_datalist, final StatsPane2 stats) {
		super();
		this.data = data;
		this.width = data[0].length * squareWidth;
		this.height = data.length * squareWidth;
		x_data = x_datalist.toArray(new Site[] {});
		y_data = y_datalist.toArray(new Site[] {});
		overall_comp=0;
		
		setPreferredSize(new Dimension(this.width, this.height));
		
		bimage = new BufferedImage(width, height, BufferedImage.TYPE_BYTE_BINARY);
		
		xmin = 0;
		ymin = 0;
		xmax = data[0].length;
		ymax = data.length;
		
		//class used for providing informative feedback when a plot is clicked
		class MousePressListener implements MouseListener {
			public void mouseClicked(MouseEvent event) {
				//get co-ordinates of mouse click
				int x = event.getX();
				int y = event.getY();
				
				//determine indices from co-ordinates
				int x_ind = 0;
				int y_ind = 0;
				if(x != 0) {
					x_ind = x / squareWidth;
				}
				if(y != 0 ) {
					y_ind = y / squareWidth;
				}
				
				//get the data from arrays, calculate percentage compatibility, output to stats pane
				if(x_ind < x_data.length && y_ind < y_data.length) {
					Site x_site = x_data[x_ind];
					Site y_site = y_data[y_ind];
					double x_percent = (x_site.getHom()/(double)(x_data.length-1))*100;
					double xp2 = Math.round(x_percent * Math.pow(10, (double) 2)) / Math.pow(10, (double) 2);
					double y_percent = (y_site.getHom()/(double)(y_data.length-1))*100;
					double yp2 = Math.round(y_percent * Math.pow(10, (double) 2)) / Math.pow(10, (double) 2);
					stats.appendln("Site "+x_data[x_ind].getPosition()+" ("+xp2+"%) cf Site "+y_data[y_ind].getPosition()+" ("+yp2+"%)");
				}
			}

			//do nothing methods
			public void mouseReleased(MouseEvent event) {}
			public void mouseEntered(MouseEvent event) {}
			public void mousePressed(MouseEvent event) {}
			public void mouseExited(MouseEvent event) {}
		}

		MouseListener listener = new MousePressListener();
		addMouseListener(listener);
		overallComp();
		drawToImage();
	}
	
	/*
		Converts plot to a graphic
	*/
	public void drawToImage() {
		Graphics2D g2i = (Graphics2D)bimage.createGraphics();
		
		//makes all of image black, so only white rectangles are drawn
		g2i.setColor(Color.black);
		g2i.fillRect(0,0,width,height);
		
		g2i.setColor(Color.white);
		
		double squareWidth = (width / data.length);
		double comp_count=0;
		
		//draws squares where sites are compatible 
		for(int i = 0; i < data.length; i++) {
			for(int j = 0; j < data[0].length; j++) {
				if(data[i][j] == true) {
					drawSquare(g2i, j, i, squareWidth);
					if(i!=j){comp_count+=1;}
				}
			}
		}
		overall_comp=comp_count/(data.length*data.length-data.length);
		g2i.dispose();
	}
	
	public void overallComp(){
		double comp_count=0;
		for(int i = 0; i < data.length; i++) {
			for(int j = 0; j < data[0].length; j++) {
				if(data[i][j] == true) {
					if(i!=j){comp_count+=1;}
				}
			}
		}
		overall_comp=comp_count/(data.length*data.length-data.length);
	}
	
	/*
		Saves image to a file.
	*/
	public void writeImageToFile(String filename, String ext) {
		try {
			ImageIO.write(bimage, ext, new File(filename + "." + ext));
		}
		catch(IOException ioe) {
			System.err.println("ERROR occurred writing image to file");
		}
	}
	public void writeImageToFile(File file, String ext) {
		try {
			ImageIO.write(bimage, ext, file);
		}
		catch(IOException ioe) {
			System.err.println("ERROR occurred writing image to file");
		}
	}
	
	/*
		Overloaded method required to draw image to panel
	*/
	public void paintComponent(Graphics g) {
		super.paintComponent(g);

		Graphics2D g2 = (Graphics2D)g;

		//draws pre-constructed image to panel
		g2.drawImage(bimage, 0, 0, null);
	}
	
	
	/*
		coordinate transformation helper methods
	*/
	public double xpixel(double xuser) {
		return (xuser - xmin) * (width) / (xmax - xmin);
	}
	public double ypixel(double yuser) {
		return (yuser - ymin) * (height) / (ymax - ymin);
	}
	
	
	/*
		Will draw a box at supplied coordinates x and y of width width
		NOTE: Will fill with whatever colour g2 has already been set to
	*/
	public void drawSquare(Graphics2D g2, double x, double y, double width) {
		Rectangle2D.Double square = new Rectangle2D.Double(xpixel(x), ypixel(y), width, width);
		g2.fill(square);
	}
	
	/*
		Returns the overall compatibility score for the alignment
		Can only be called once are grid has been 'drawn' up
	*/
	public double getOverallComp()
	{
		return overall_comp;
	}
}