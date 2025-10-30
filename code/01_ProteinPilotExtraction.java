import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
//import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
//import java.util.Map.Entry;
import java.util.Set;

import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.DateUtil;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;


/**
 * Hello world!
 *
 */
public class ProteinPilotExtraction 
{
	
	private static String correctingMods(String mods, int pepstart, String pep, String accession) {
		String newmodslist = "";
		String[] protmodlist = mods.split(";");
		for(int i = 0; i < protmodlist.length; ++i) {
			String mod = protmodlist[i].trim();
			int idx = mod.indexOf("@");
			String mod_name = "";
			if(idx>-1) {
				mod_name= mod.substring(0, idx);
				int idx2 = mod_name.indexOf(")");
				if(idx2>-1) {
					mod_name = mod_name.substring(0, idx2+1);
				}
				String siteidx = mod.substring(idx+1);
				int siteidx_i = -1;
				if(siteidx.startsWith("N-term")) {
					siteidx_i = 0;
				} else if(siteidx.startsWith("C-term")){
					siteidx_i = pep.length();
				} else {
					siteidx_i = Integer.parseInt(siteidx);
				}
				siteidx_i = siteidx_i + pepstart;
				
				if(mod_name.equals("Phospho(S)") || mod_name.equals("Phospho(T)") || mod_name.equals("Phospho(Y)") || 
						/*mod_name.equals("Acetyl") || */mod_name.equals("Acetyl(K)") || mod_name.equals("GG(K)") || mod_name.equals("GG(S)") ||
						mod_name.equals("Methyl(K)") || mod_name.equals("Methyl(R)") || mod_name.equals("Oxidation(M)") || 
						mod_name.equals("Deamidated(N)") || mod_name.equals("Deamidated(Q)") || mod_name.equals("Oxidation(P)") ||
						mod_name.equals("Deamidated(R)") || mod_name.equals("Citrullination(R)")) {
					if(mod_name.equals("Deamidated(R)")) {
						mod_name = "Citrullination(R)";
					}
					if(accession.equals("3R")) {
						if(siteidx_i>20) {
							siteidx_i = siteidx_i + 31;
						}
						siteidx_i = siteidx_i + 254;
					} else if(accession.equals("0N")) {
						if(siteidx_i>20) {
							siteidx_i = siteidx_i + 58;
						}
						siteidx_i = siteidx_i + 24;
					} else if(accession.equals("1N")) {
						if(siteidx_i>49) {
							siteidx_i = siteidx_i + 29;
						}
						siteidx_i = siteidx_i + 24;
					} else if(accession.equals("epsilon")) {
						siteidx_i = siteidx_i + 376;
					} else if(accession.equals("Ex3")) {
						if(siteidx_i>7) {
							siteidx_i = siteidx_i + 14;
						}
						siteidx_i = siteidx_i + 32;
					} else if(accession.equals("Ex5")) {
						if(siteidx_i>7) {
							siteidx_i = siteidx_i + 28;
						}
						siteidx_i = siteidx_i + 96;
					}
					String tmp = mod_name + "@" + siteidx_i;
					if(pepstart>-1) {
						if(!newmodslist.equals("")) {
							newmodslist = newmodslist + "; ";
						}
						newmodslist = newmodslist + tmp;
					}
				}
			}
		}
		return newmodslist;
	}
	
	private static String correctingCleavages(String cleavages, int pepstart, String pep, String accession) {
		String newcleavagelist = "";
		String[] cleavagelist = cleavages.split(";");
		for(int i = 0; i < cleavagelist.length; ++i) {
			String mod = cleavagelist[i].trim();
			if(mod.startsWith("cleaved")) {
				int idx = mod.indexOf("@");
//				String cleave_name = "";
				if(idx>-1) {
//					cleave_name= mod.substring(0, idx);
					String siteidx = mod.substring(idx+1);
					int siteidx_i = -1;
					String aminoacid = "";
					if(siteidx.startsWith("N-term")) {
						siteidx_i = 1;
						aminoacid = pep.substring(0,1);
					} else if(siteidx.startsWith("C-term")){
						siteidx_i = pep.length();
						aminoacid = pep.substring(pep.length()-1);
					}
					siteidx_i = siteidx_i + pepstart;
					
					if(siteidx.equals("N-term") || siteidx.equals("C-term")) {
						if(accession.equals("3R")) {
							if(siteidx_i>20) {
								siteidx_i = siteidx_i + 31;
							}
							siteidx_i = siteidx_i + 254;
						} else if(accession.equals("0N")) {
							if(siteidx_i>20) {
								siteidx_i = siteidx_i + 58;
							}
							siteidx_i = siteidx_i + 24;
						} else if(accession.equals("1N")) {
							if(siteidx_i>49) {
								siteidx_i = siteidx_i + 29;
							}
							siteidx_i = siteidx_i + 24;
						} else if(accession.equals("epsilon")) {
							siteidx_i = siteidx_i + 376;
						} else if(accession.equals("Ex3")) {
							if(siteidx_i>7) {
								siteidx_i = siteidx_i + 14;
							}
							siteidx_i = siteidx_i + 32;
						} else if(accession.equals("Ex5")) {
							if(siteidx_i>7) {
								siteidx_i = siteidx_i + 28;
							}
							siteidx_i = siteidx_i + 96;
						}
						String tmp = siteidx + "(" + aminoacid + ")@" + siteidx_i;
						if(pepstart>-1) {
							if(!newcleavagelist.equals("")) {
								newcleavagelist = newcleavagelist + "; ";
							}
							newcleavagelist = newcleavagelist + tmp;
						}
					}
				}
			}
		}
		return newcleavagelist;
	}
	
	private static boolean checkmods(String ProtModifications) {
		boolean result = false;
		
		Set<String> protmodtable_mods = new HashSet<String>();
		
		String[] protmodlist = ProtModifications.split(";");
		for(int i = 0; i < protmodlist.length; ++i) {
			String mod = protmodlist[i].trim();
			int idx = mod.indexOf("@");
			String mod_name = "";
			if(idx>-1) {
				mod_name= mod.substring(0, idx);
				int idx2 = mod_name.indexOf(")");
				if(idx2>-1) {
					mod_name = mod_name.substring(0, idx2+1);
				}
				if(mod_name.equals("Deamidated(R)")) {
					mod_name = "Citrullination(R)";
				}
				protmodtable_mods.add(mod_name);
			}
		}
		int count = 0;
		int count_interest = 0;
		if(protmodtable_mods.contains("Phospho(S)") || protmodtable_mods.contains("Phospho(T)") || protmodtable_mods.contains("Phospho(Y)")) {
			count += 1;
			count_interest += 1;
		}
		if (protmodtable_mods.contains("GG(K)") || protmodtable_mods.contains("GG(S)")) {
			count += 1;
			count_interest += 1;
		}
		if (/*protmodtable_mods.contains("Acetyl") || */protmodtable_mods.contains("Acetyl(K)")) {
			count += 1;
			count_interest += 1;
		}
		if (protmodtable_mods.contains("Methyl(K)") | protmodtable_mods.contains("Methyl(R)")) {
			count += 1;
			count_interest += 1;
		}
		if (protmodtable_mods.contains("Citrullination(R)") | protmodtable_mods.contains("Deamidated(R)")) {
			count += 1;
			count_interest += 1;
		}
		if (protmodtable_mods.contains("Deamidated(N)")) {
			count += 1;
			count_interest += 1;
		}
		if (protmodtable_mods.contains("Oxidation(P)")) {
			count += 1;
			count_interest += 1;
		}
		if (protmodtable_mods.contains("Oxidation(M)")) {
			count += 1;
		}
		if (protmodtable_mods.contains("Deamidated(Q)")) {
			count += 1;
		}
		
		if(count==protmodtable_mods.size() && count_interest>0) {
			result = true;
		} else {
			result = false;
		}
		
		return result;
	}
	
	private static boolean checkcleavages(String cleavages) {
		boolean result = false;
		
		String[] cleavagelist = cleavages.split(";");
		for(int i = 0; i < cleavagelist.length; ++i) {
			String mod = cleavagelist[i].trim();
			int idx = mod.indexOf("@");
			String mod_name = "";
			if(idx>-1) {
				mod_name= mod.substring(0, idx);
				if(mod_name.startsWith("cleaved")) {
					result = true;
				}
			}
		}
		
		return result;
	}
	
	
	private static String addPTMsequences(String sequence, String modifications, String ProtModifications) {
		List<String> protmodtable_mods = new ArrayList<String>();
		List<String> protmodtable_name = new ArrayList<String>();
		List<Integer> protmodtable_idxs = new ArrayList<Integer>();
		String modsequence = "";
//		System.out.println(sequence + "\t" + modifications + "\t" + ProtModifications);
		if(!sequence.equals("") && !modifications.equals("") && !ProtModifications.equals("")) {
			modsequence = sequence;
			String[] modlist = modifications.split(";");
			String[] protmodlist = ProtModifications.split(";");
			for(int j = protmodlist.length-1; j >= 0; --j) {
				String mod = protmodlist[j].trim();
				if(!mod.contains("-term")) {
					int idx = mod.indexOf("(");
					String mod_name = "";
					if(idx>-1) {
						mod_name= mod.substring(0, idx);
						int idx2 = mod.indexOf(")");
						String mod_site = mod.substring(idx+1, idx2);
						int idx3 = mod.indexOf("@");
						int mod_index = Integer.parseInt(mod.substring(idx3+1)) - 1;
						protmodtable_mods.add(mod_name + "(" + mod_site + ")");
						protmodtable_name.add(mod_name);
						protmodtable_idxs.add(mod_index);
					}
					
				} else {
					int idx = mod.indexOf("@");
					String mod_name = mod.substring(0, idx);
					String mod_site = mod.substring(idx+1);
					if(mod_site.equals("N-term")) {
						protmodtable_mods.add(mod_name + "(" + mod_site + ")");
						protmodtable_name.add(mod_name);
						protmodtable_idxs.add(-1);
					} else if(mod_site.equals("C-term")) {
						protmodtable_mods.add(mod_name + "(" + mod_site + ")");
						protmodtable_name.add(mod_name);
						protmodtable_idxs.add(sequence.length()+1);
					}
				}
			}			
			
			for(int i = modlist.length-1; i >= 0; --i) {
				String mod = modlist[i].trim();
				if(!mod.contains("-term")) {
					int idx = mod.indexOf("(");
					String mod_name = mod.substring(0, idx);
					int idx2 = mod.indexOf(")");
					String mod_site = mod.substring(idx+1, idx2);
					int idx3 = mod.indexOf("@");
					int mod_index = Integer.parseInt(mod.substring(idx3+1)) - 1;
					boolean found = false;
					int protmodidx = 0;
					while(found!=true && protmodidx<protmodtable_mods.size()) {
						if(protmodtable_name.get(protmodidx).equals(mod_name)) {
							if(protmodtable_mods.get(protmodidx).equals(mod_name + "(" + mod_site + ")")) {
								found = true;
							} else if(mod_site.equals("N-term") || mod_site.equals("C-term")){
								found = true;
							}
						}
						protmodidx += 1;
					}
					if(found==true) {
						protmodidx = protmodidx - 1;
						String pepsite = "" + modsequence.charAt(mod_index);
						if(mod_site.equals(pepsite)) {
							modsequence = modsequence.substring(0, mod_index + 1) + "(" + mod_name + ")" + modsequence.substring(mod_index + 1);
//							String protmod_mod = protmodtable_mods.remove(protmodidx);
//							String protmod_nam = protmodtable_name.remove(protmodidx);
//							Integer protmod_idx = protmodtable_idxs.remove(protmodidx);
						}
					}
				} else {
					int idx = mod.indexOf("@");
					String mod_name = mod.substring(0, idx);
					String mod_site = mod.substring(idx+1);
					boolean found = false;
					int protmodidx = 0;
					while(found!=true && protmodidx<protmodtable_mods.size()) {
						if(protmodtable_name.get(protmodidx).equals(mod_name)) {
							if(protmodtable_mods.get(protmodidx).equals(mod_name + "(" + mod_site + ")")) {
								found = true;
							} else if(mod_site.equals("N-term") || mod_site.equals("C-term")){
								found = true;
							}
						}
						protmodidx += 1;
					}
					if(found==true) {
						protmodidx = protmodidx - 1;
						if(mod_site.equals("N-term")) {
							modsequence = "(" + mod_name + ")" + modsequence;
//							String protmod_mod = protmodtable_mods.remove(protmodidx);
//							String protmod_nam = protmodtable_name.remove(protmodidx);
//							Integer protmod_idx = protmodtable_idxs.remove(protmodidx);
						} else if(mod_site.equals("C-term")) {
							modsequence = modsequence + "(" + mod_name + ")";
//							String protmod_mod = protmodtable_mods.remove(protmodidx);
//							String protmod_nam = protmodtable_name.remove(protmodidx);
//							Integer protmod_idx = protmodtable_idxs.remove(protmodidx);
						}
					}
				}
			}
		} else if(!sequence.equals("") && modifications.equals("")) {
			modsequence = sequence;
		}
		return modsequence;
	}
	
//	private static String idincluded(String testing, String[] inputids, List<String> headers_ordered) {
//		String output = "";
//		int idx = 0;
//		boolean found = false;
//		String idtest = "";
//		while(idx<inputids.length && found==false) {
//			idtest = inputids[idx];
//			found = testing.contains(idtest);
//			idx += 1;
//		}
//		if(found==true) {
//			found = false;
//			idx = 0;
//			String header = "";
//			while(idx<headers_ordered.size() && found==false) {
//				header = headers_ordered.get(idx);
//				found = header.contains(idtest);
//				idx += 1;
//			}
//			
//			if(found==true) {
//				output = header;
//			}
//		}
//		return output;
//	}
	
	
    @SuppressWarnings("resource")
	public static void main( String[] args )
    {
    	String taufile = args[0]; //"D:/Data/Simon_HMW_paper/Addition/uniprot_tau_isoforms2N4R-0N-1N-3R.fasta"; // 
    	List<String> headers_ordered = new ArrayList<String>();
    	Map<String, String> taumap = new HashMap<String,String>();
    	
    	BufferedReader reader;
    	try {
    		reader = new BufferedReader(new FileReader(new File(taufile)));
    		String line = "";
    		String header = "";
    		String sequence = "";
    		while((line=reader.readLine())!=null) {
    			if(line.startsWith(">")) {
    				if(!header.equals("")) {
    					headers_ordered.add(header);
    					taumap.put(header, sequence);
    					header = "";
    					sequence = "";
    				}
    				header = line.split("\\|")[1];
    			} else {
    				sequence = sequence + line.trim();
    			}
    		}
    		headers_ordered.add(header);
			taumap.put(header, sequence);
    		reader.close();
    	} catch (IOException e) {
    		// TODO auto-generated catch block
    		e.printStackTrace();
    	}
    	
//    	System.out.println(headers_ordered);
//    	System.out.println(taumap);
//    	System.out.println(taumap.get(headers_ordered.get(0)) + "\n" + "AKTDHGAEIVYK");
//    	System.out.println(taumap.get(headers_ordered.get(0)).contains("AKTDHGAEIVYK"));
//    	System.out.println(taumap.get(headers_ordered.get(0)).indexOf("AKTDHGAEIVYK"));
    	
    	
        String filetablelocation = args[1]; //"D:/Data/Simon_HMW_paper/Addition/00_FDR_files_HMW_CTR.txt"; // 
    	String path = args[2]; //"D:/Data/Simon_HMW_paper/Addition"; // 
    	String outputfile = args[3]; //"D:/Data/Simon_HMW_paper/Addition/PP_cutoff_0-05FDR_HMW_CTR.csv"; // 
    	double threshold = Double.parseDouble(args[4]); //1; //0.05; 
    	String id = args[5]; //"P10636"; //
//    	String ids_string = args[5]; //"P10636"; //
//    	String[] ids = ids_string.split(",");
    	int fdrtype = (int)Double.parseDouble(args[6]); // 17 == global FDR; 18 == local FDR
    	String sep = ",";
    	List<String> files = new ArrayList<String>();
    	
    	try {
			reader = new BufferedReader(new FileReader(new File(filetablelocation)));
			String line = "";
			while((line=reader.readLine())!=null) {
				files.add(line);
			}
			reader.close();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
    	
    	Writer writer;
		try {
			writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outputfile)));
//			writer.write("Groupfile" + sep + "FDR" + sep + "FDR Threshold" + sep + "Conf Threshold\n");
			writer.write("Sample" + sep + "Accessions" + sep + "Conf" + sep + "Sequence" + sep + "Modificaitons" + sep + "CorrectedModifications" + sep + "ProteinModifications" + sep + "ModifiedSequence" + sep + "Cleavages" + sep + "CleavageSites" + sep + "dMass" + sep + "Obs MW" + sep + "Obs m/z" + sep + "Theor MW" + sep + "Theor m/z" + sep + "Theor z" + "\n");
			
			for(String file : files) {
		    	
		    	String fileLocation = path + File.separatorChar + file;
		    	file = file.replace("__FDR.xlsx", "_PeptideSummary.txt");
		    	String[] filesplit = file.split("/");
		    	String newfile = filesplit[filesplit.length-1];
		    	String[] newfilesplit = newfile.split("_");
		    	String newfilename = newfilesplit[0] + "_" + newfilesplit[1] + "_" + newfilesplit[2];// + "_" + newfilesplit[3];
//		    	String newfilename = newfilesplit[2] + "_" + newfilesplit[4];
				FileInputStream xlsxfile = new FileInputStream(new File(fileLocation));
				Workbook workbook = new XSSFWorkbook(xlsxfile);
				int idx = 0;
				
				for(Sheet sheets : workbook) {
					if(sheets.getSheetName().equals("Protein Level Data")) {
						break;
					}
					idx += 1;
				}
				
				Sheet sheet = workbook.getSheetAt(idx);
				//System.out.println(sheet.getSheetName());
				List<Integer> myarray = new ArrayList<Integer>();
				myarray.add(16);
				myarray.add(17);
				myarray.add(18);
				boolean Global1found = false;
				boolean Global5found = false;
				boolean Local1found = false;
				boolean Local5found = false;
//				double FDR = 0.0;
				int Global1number = 0;
				int Global5number = 0;
				int Local1number = 0;
				int Local5number = 0;
				for (Row row : sheet) {
					for(int i : myarray) {
						Cell cell = row.getCell(i);
						if(cell!=null) {
					    	switch (cell.getCellType()) {
					            case STRING: break;
					            case NUMERIC: if (DateUtil.isCellDateFormatted(cell)) {
					            					break;
									            } else {
									            	double value = cell.getNumericCellValue();
									            	if(i==18) {
//									            		System.out.println("FDR=" + value + "<" + threshold + "==" + (value > threshold));
									            		if(value>=0.01 && Local1found!=true) {
									            			Local1found = true;
									            			if(Local5found==true && Global1found==true && Global5found==true) {
									            				break;
									            			}
									            		}
									            		if(value>=0.05 && Local5found!=true) {
									            			Local5found = true;
									            			if(Local1found==true && Global1found==true && Global5found==true) {
									            				break;
									            			}
									            		}
									            	} else if(i==17) {
									            		if(value>=0.01 && Global1found!=true) {
									            			Global1found = true;
									            			if(Local1found==true && Local5found==true && Global5found==true) {
									            				break;
									            			}
									            		}
									            		if(value>=0.05 && Global5found!=true) {
									            			Global5found = true;
									            			if(Local1found==true && Local5found==true && Global1found==true) {
									            				break;
									            			}
									            		}
									            	
									            	} else if (i==16) {
//									            		System.out.println("Conf=" + value);
									            		if(Local1found!=true) {
									            			Local1number = (int)value;
									            		}
									            		if(Local5found!=true) {
									            			Local5number = (int)value;
									            		}
									            		if(Global1found!=true) {
									            			Global1number = (int)value;
									            		}
									            		if(Global5found!=true) {
									            			Global5number = (int)value;
									            		}
									            	}
									            } break;
					            case BOOLEAN: break;
					            case FORMULA: break;
					            default: break;
					    	}
						} else {
							break;
						}
				    	if(Local1found==true && Local5found==true && Global1found==true && Global5found==true) {
				    		break;
				    	}
				    }
					if(Local1found==true && Local5found==true && Global1found==true && Global5found==true) {
						break;
					}
				}
				
				
				
				idx = 0;
				
				for(Sheet sheets : workbook) {
					if(sheets.getSheetName().equals("Protein Summary")) {
						break;
					}
					idx += 1;
				}
				
				sheet = workbook.getSheetAt(idx);
				//System.out.println(sheet.getSheetName());
				myarray = new ArrayList<Integer>();
				myarray.add(0);
				myarray.add(6);
				boolean Taufound = false;
				boolean Syafound = false;
				boolean Tdpfound = false;
				boolean Appfound = false;
				boolean Padi2found = false;
				boolean Gfapfound = false;
//				double FDR = 0.0;
				int Taurank = 0;
				int Syarank = 0;
				int Tdprank = 0;
				int Apprank = 0;
				int Padi2rank = 0;
				int Gfaprank = 0;
				for (Row row : sheet) {
					for(int i : myarray) {
						Cell cell = row.getCell(i);
						if(cell!=null) {
					    	switch (cell.getCellType()) {
					            case STRING: String value = cell.getStringCellValue();
								            	if(i==6) {
				//				            		System.out.println("FDR=" + value + "<" + threshold + "==" + (value > threshold));
								            		if(value.contains("P10636") && Taufound!=true) {
								            			Taufound = true;
								            			if(Syafound==true && Tdpfound==true && Appfound==true && Padi2found==true && Gfapfound==true) {
								            				break;
								            			}
				//				            		} else {
				////				            			FDR = value;
								            		}
								            		if(value.contains("P37840") && Syafound!=true) {
								            			Syafound = true;
								            			if(Taufound==true && Tdpfound==true && Appfound==true && Padi2found==true && Gfapfound==true) {
								            				break;
								            			}
								            		}
								            		if(value.contains("Q13148") && Tdpfound!=true) {
								            			Tdpfound = true;
								            			if(Taufound==true && Syafound==true && Appfound==true && Padi2found==true && Gfapfound==true) {
								            				break;
								            			}
								            		}
								            		if(value.contains("P05067") && Appfound!=true) {
								            			Appfound = true;
								            			if(Taufound==true && Syafound==true && Tdpfound==true && Padi2found==true && Gfapfound==true) {
								            				break;
								            			}
								            		}
								            		if(value.contains("Q9Y2J8") && Padi2found!=true) {
								            			Padi2found = true;
								            			if(Taufound==true && Syafound==true && Tdpfound==true && Appfound==true && Gfapfound==true) {
								            				break;
								            			}
								            		}
								            		if(value.contains("P14136") && Gfapfound!=true) {
								            			Gfapfound = true;
								            			if(Taufound==true && Syafound==true && Tdpfound==true && Appfound==true && Padi2found==true) {
								            				break;
								            			}
								            		}
								            	} break;
					            case NUMERIC: if (DateUtil.isCellDateFormatted(cell)) {
					            					break;
									            } else {
									            	int valueint = (int)cell.getNumericCellValue();
									            	if (i==0) {
//									            		System.out.println("Conf=" + value);
									            		if(Taufound!=true) {
									            			Taurank = valueint;
									            		}
									            		if(Syafound!=true) {
									            			Syarank = valueint;
									            		}
									            		if(Tdpfound!=true) {
									            			Tdprank = valueint;
									            		}
									            		if(Appfound!=true) {
									            			Apprank = valueint;
									            		}
									            		if(Padi2found!=true) {
									            			Padi2rank = valueint;
									            		}
									            		if(Gfapfound!=true) {
									            			Gfaprank = valueint;
									            		}
									            	}
									            } break;
					            case BOOLEAN: break;
					            case FORMULA: break;
					            default: break;
					    	}
						} else {
							break;
						}
				    	if(Taufound==true && Syafound==true && Tdpfound==true && Appfound==true && Padi2found==true && Gfapfound==true) {
				    		break;
				    	}
				    }
					if(Taufound==true && Syafound==true && Tdpfound==true && Appfound==true && Padi2found==true && Gfapfound==true) {
						break;
					}
				}
				
				System.out.println(file + "\tglob1%=" + Global1number + "\tglob5%=" + Global5number + "\tloc1%=" + Local1number + "\tloc5%=" + Local5number + "\tTau=" + Taurank + "\tSyn=" + Syarank + "\tTdp=" + Tdprank + "\tApp=" + Apprank + "\tPadi2=" + Padi2rank + "\tGfap=" + Gfaprank);
				
				idx = 0;
				
				for(Sheet sheets : workbook) {
					if(sheets.getSheetName().equals("Peptide Summary")) {
						break;
					}
					idx += 1;
				}
				
				sheet = workbook.getSheetAt(idx);
				//System.out.println(sheet.getSheetName());

//				List<Map<String, String>> list = new ArrayList<Map<String, String>>();
				 // Get the maximum number of rows
				int rownum = sheet.getPhysicalNumberOfRows();
				 // Get the second line
				Row row_tmp = sheet.getRow(0);
				 // Get the maximum number of columns
				int colnum = row_tmp.getPhysicalNumberOfCells();
//				String columns[] = {"N","Unused","Total","%Cov","%Cov(50)","%Cov(95)","Accessions","Names","Used","Annotation","Contrib","Conf","Sequence","Modifications","ProteinModifications","Cleavages","dMass","Obs MW","Obs m/z","Theor MW","Theor m/z","Theor z","Sc","Spectrum","Acq Time","Intensity (Peptide)","PrecursorIntensityAcquisition","Apex Time (Peptide)","Elution Peak Width (Peptide)","MS2Counts"};
				String cellData = null;
				
				File tmpfile = new File(path + File.separatorChar + file);
				if(!tmpfile.exists()) {
					Writer writer_tmp = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(path + File.separatorChar + file)));
					
					for (int i = 0; i < rownum; i++) {
//						Map<String, String> map = new LinkedHashMap<String, String>();
						row_tmp = sheet.getRow(i);
						if (row_tmp != null) {
							for (int j = 0; j < colnum; j++) {
								if(j!=0) {
									writer_tmp.write("\t");
								}
								cellData = (String) getCellFormatValue(row_tmp.getCell(j));
//								map.put(columns[j], cellData);
								if(cellData==null) {
									cellData="";
								}
								writer_tmp.write(cellData);
							}
							writer_tmp.write("\r\n");
						} else {
							break;
						}
//						list.add(map);
					}
	
					writer_tmp.close();
				}
				// Traverse the parsed list
//				StringBuffer sb = new StringBuffer();
//				for (int i = 0; i < list.size(); i++) {
//					boolean start = true;
//					for (Entry<String, String> entry : list.get(i).entrySet()) {
//						if(start==false) {
//							sb.append("\t");
//						}
//						String value = entry.getValue();
//						sb.append(value);
//						start=false;
//					}
//					sb.append("\r\n");
//				}
//				WriteToFile(sb.toString(), path + File.separatorChar + file);
				
				idx = 0;
				for(Sheet sheets : workbook) {
					if(sheets.getSheetName().equals("Distinct Peptide Level Data")) {
						break;
					}
					idx += 1;
				}
				
				sheet = workbook.getSheetAt(idx);
				//System.out.println(sheet.getSheetName());
				myarray = new ArrayList<Integer>();
				myarray.add(fdrtype); // 17 == global FDR; 18 == local FDR
				myarray.add(11);
				boolean found = false;
//				double FDR = 0.0;
				double conf = 0.0;
				for (Row row : sheet) {
					for(int i : myarray) {
						Cell cell = row.getCell(i);
						if(cell!=null) {
					    	switch (cell.getCellType()) {
					            case STRING: break;
					            case NUMERIC: if (DateUtil.isCellDateFormatted(cell)) {
					            					break;
									            } else {
									            	double value = cell.getNumericCellValue();
									            	if(i==fdrtype) {
//									            		System.out.println("FDR=" + value + "<" + threshold + "==" + (value > threshold));
									            		if(value > threshold) {
									            			found = true;
									            			break;
//									            		} else {
////									            			FDR = value;
									            		}
									            	} else if (i==11) {
//									            		System.out.println("Conf=" + value);
									            		conf = value;
									            	}
									            } break;
					            case BOOLEAN: break;
					            case FORMULA: break;
					            default: break;
					    	}
						} else {
							break;
						}
				    	if(found==true) {
				    		break;
				    	}
				    }
					if(found==true) {
						break;
					}
				}
				
				xlsxfile.close();
				
				BufferedReader pepsumreader = new BufferedReader(new FileReader(new File(path + File.separatorChar + file)));
				String pepsumline = "";
				while((pepsumline=pepsumreader.readLine())!=null) {
					if(!pepsumline.startsWith("N")) {
						String[] pepsumlinesplit = pepsumline.split("\t");
						if(pepsumlinesplit.length>0) {
//							String id = idincluded(pepsumlinesplit[6],ids,headers_ordered);
//							if(!id.equals("")) {
							if(pepsumlinesplit[6].contains(id)) {
								double pepsumlineconf = Double.parseDouble(pepsumlinesplit[11]);
//								System.out.println(pepsumlineconf + ">=" + (conf*100.0) + "==" + (pepsumlineconf >= (conf*100.0)));
								if(pepsumlineconf >= (conf*100.0)) {
									String accession = "";
									int index = -1;
//									if(headers_ordered.size()>4) {
//										if(taumap.get(id).contains(pepsumlinesplit[12])) {
//											accession = id;
//											index = taumap.get(id).indexOf(pepsumlinesplit[12]);
//										}
//									} else 
										if(headers_ordered.size()==4) {
										if(taumap.get(headers_ordered.get(0)).contains(pepsumlinesplit[12])) {
											accession = "2N4R";
											index = taumap.get(headers_ordered.get(0)).indexOf(pepsumlinesplit[12]);
										} else if(taumap.get(headers_ordered.get(1)).contains(pepsumlinesplit[12])) {
											accession = "0N";
											index = taumap.get(headers_ordered.get(1)).indexOf(pepsumlinesplit[12]);
										} else if(taumap.get(headers_ordered.get(2)).contains(pepsumlinesplit[12])) {
											accession = "1N";
											index = taumap.get(headers_ordered.get(2)).indexOf(pepsumlinesplit[12]);
										} else if(taumap.get(headers_ordered.get(3)).contains(pepsumlinesplit[12])) {
											accession = "3R";
											index = taumap.get(headers_ordered.get(3)).indexOf(pepsumlinesplit[12]);
										}
									} else if(headers_ordered.size()==2) {
										if(taumap.get(headers_ordered.get(0)).contains(pepsumlinesplit[12])) {
											accession = "alpha";
											index = taumap.get(headers_ordered.get(0)).indexOf(pepsumlinesplit[12]);
										} else if(taumap.get(headers_ordered.get(1)).contains(pepsumlinesplit[12])) {
											accession = "epsilon";
											index = taumap.get(headers_ordered.get(1)).indexOf(pepsumlinesplit[12]);
										}
									} else if(headers_ordered.size()==3) {
										if(taumap.get(headers_ordered.get(0)).contains(pepsumlinesplit[12])) {
											accession = "Can";
											index = taumap.get(headers_ordered.get(0)).indexOf(pepsumlinesplit[12]);
										} else if(taumap.get(headers_ordered.get(1)).contains(pepsumlinesplit[12])) {
											accession = "Ex5";
											index = taumap.get(headers_ordered.get(1)).indexOf(pepsumlinesplit[12]);
										} else if(taumap.get(headers_ordered.get(2)).contains(pepsumlinesplit[12])) {
											accession = "Ex3";
											index = taumap.get(headers_ordered.get(2)).indexOf(pepsumlinesplit[12]);
										}
									} else if(headers_ordered.size()==1) {
										if(taumap.get(headers_ordered.get(0)).contains(pepsumlinesplit[12])) {
											accession = "PADI2";
											index = taumap.get(headers_ordered.get(0)).indexOf(pepsumlinesplit[12]);
										}
									}
	//								if(newfilename.equals("PSP18_BA39-Cohort1_INSOLUBLE")) {
	//									if(pepsumlinesplit[12].equals("AGLKAEEAGIGDTPSLEDEAAGHVTQAR")) {
	//										if(pepsumlinesplit[14].equals("Phospho(T)@53; Oxidation(P)@54; Phospho(S)@55")) {
	//											boolean test = true;
	//											System.out.println(test);
	//										}
	//									}
	//								}
//									System.out.println(pepsumlinesplit[12]);
//									if(pepsumlinesplit[12].startsWith("VIPIQAHQIV")) {
//										System.out.println("checkmods: " + checkmods(pepsumlinesplit[14]));
//										System.out.println("checkcleavages: " + checkcleavages(pepsumlinesplit[15]));
//										System.out.println("checkpepmods: " + (pepsumlinesplit[14].equals("")==true && checkmods(pepsumlinesplit[13])==true));
//									}
									if((checkmods(pepsumlinesplit[14])==true || checkcleavages(pepsumlinesplit[15])==true || (pepsumlinesplit[14].equals("")==true && checkmods(pepsumlinesplit[13])==true)) || (pepsumlinesplit[14].equals("")==true && pepsumlinesplit[13].equals("")==true) && !accession.equals("")) {
										writer.write(newfilename + sep + accession + sep + pepsumlinesplit[11] + sep + pepsumlinesplit[12] + sep + pepsumlinesplit[13] + sep + correctingMods(pepsumlinesplit[13], index, pepsumlinesplit[12], accession) + sep);
										writer.write(pepsumlinesplit[14] + sep + addPTMsequences(pepsumlinesplit[12], pepsumlinesplit[13], pepsumlinesplit[14]) + sep);
										writer.write(pepsumlinesplit[15] + sep + correctingCleavages(pepsumlinesplit[15], index, pepsumlinesplit[12], accession) + sep + pepsumlinesplit[16] + sep + pepsumlinesplit[17] + sep);
										writer.write(pepsumlinesplit[18] + sep + pepsumlinesplit[19] + sep + pepsumlinesplit[20] + sep + pepsumlinesplit[21] + "\n");
									}
								}
							}
						}
					}
				}
				
				pepsumreader.close();
//				writer.write(newfile + sep + FDR + sep + threshold + sep + conf + "\n");
	    	}
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    }
    
    public static void WriteToFile(String str, String filePath) throws IOException {
		BufferedWriter bw = null;
		try {
			 FileOutputStream out = new FileOutputStream(filePath, true); // true, means: file appended content, not regenerated, default is false
			bw = new BufferedWriter(new OutputStreamWriter(out));
			 bw.write(str += "\r\n");// Line break
			bw.flush();
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			bw.close();
		}
    }

	public static Object getCellFormatValue(Cell cell) {
		Object cellValue = null;
		if (cell != null) {
			 // Determine the cell type
			switch (cell.getCellType()) {
			case NUMERIC: {
				String val = String.valueOf(cell.getNumericCellValue());
				int idx = val.indexOf("E");
				if(idx>-1) {
					Double exponent = Double.parseDouble(val.substring(idx+1));
					if(exponent>=7.0) {
						String tmpval = String.valueOf(String.format("%.6E",cell.getNumericCellValue()));
						int tmpidx = tmpval.indexOf("0E");
						int count = 0;
						while(tmpidx>-1) {
							count += 1;
							tmpval = tmpval.replace("0E","E");
							tmpidx = tmpval.indexOf("0E");
						}
						if(count==0) {
							cellValue = String.valueOf(String.format("%.6E",cell.getNumericCellValue()));
						} else if(count==1) {
							cellValue = String.valueOf(String.format("%.5E",cell.getNumericCellValue()));
						} else if(count==2) {
							cellValue = String.valueOf(String.format("%.4E",cell.getNumericCellValue()));
						} else if(count==3) {
							cellValue = String.valueOf(String.format("%.3E",cell.getNumericCellValue()));
						}
					} else if(exponent<=-5.0) {	
						String tmpval = String.valueOf(String.format("%.14E",cell.getNumericCellValue()));
						int tmpidx = tmpval.indexOf("0E");
						int count = 0;
						while(tmpidx>-1) {
							count += 1;
							tmpval = tmpval.replace("0E","E");
							tmpidx = tmpval.indexOf("0E");
						}
						if(count==0) {
							cellValue = String.valueOf(String.format("%.14E",cell.getNumericCellValue()));
						} else if(count==1) {
							cellValue = String.valueOf(String.format("%.13E",cell.getNumericCellValue()));
						} else if(count==2) {
							cellValue = String.valueOf(String.format("%.12E",cell.getNumericCellValue()));
						} else if(count==3) {
							cellValue = String.valueOf(String.format("%.11E",cell.getNumericCellValue()));
						} else if(count==4) {
							cellValue = String.valueOf(String.format("%.10E",cell.getNumericCellValue()));
						}
					} else {
						DecimalFormat df = new DecimalFormat("0");
						df.setMaximumFractionDigits(340);
						cellValue = String.valueOf(df.format(cell.getNumericCellValue()));
					}
				} else {
					DecimalFormat df = new DecimalFormat("0");
					df.setMaximumFractionDigits(340);
					cellValue = String.valueOf(df.format(cell.getNumericCellValue()));
				}
				
				break;
			}
			case FORMULA: {
				 // Determine whether the cell is in date format
				if (DateUtil.isCellDateFormatted(cell)) {
					 // Convert to date format YYYY-mm-dd
					cellValue = cell.getDateCellValue();
				} else {
					 // number
					cellValue = String.valueOf(cell.getNumericCellValue());
				}
				break;
			}
			case STRING: {
				cellValue = cell.getRichStringCellValue().getString();
				break;
			}
			case BOOLEAN: {
				cellValue = cell.getBooleanCellValue();
			}
			default:
				cellValue = "";
			}
		} else {
			cellValue = "";
		}
		return cellValue;
	}    
}
