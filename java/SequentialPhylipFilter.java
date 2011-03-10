import java.io.File;
import javax.swing.*;
import javax.swing.filechooser.*;

public class SequentialPhylipFilter extends FileFilter {
	public boolean accept(File f) {
		return true;
	}
	public String getDescription() {
		return "Sequential Phylip Files";
	}
}