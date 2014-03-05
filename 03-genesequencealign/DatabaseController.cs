using System;
using System.Collections.Generic;
using System.Text;
using System.Data.OleDb;
using System.IO;

namespace GeneticsLab
{
    class DatabaseController
    {
        private OleDbConnection m_accessConn = null;	// connection to the database

        public DatabaseController()
        {
        }

        public void EstablishConnection(string dbaseFileName)
        {
            string strAccessConn = "Provider=Microsoft.Jet.OLEDB.4.0;Data Source=" + dbaseFileName;
            m_accessConn = new OleDbConnection(strAccessConn);
        }

        public GeneSequence ReadGeneSequence(int id)
        {
            int problemCount;
            string name = "[error]", sequence = "[error]"; // default values in case of a exception

            try
            {
                OleDbCommand countSequencesCommand = new OleDbCommand("SELECT MAX(ID) FROM DNA", m_accessConn);
                //OleDbCommand selectCommand = new OleDbCommand("SELECT * FROM DNA WHERE ID = " + id, m_accessConn);
                OleDbCommand selectCommand = new OleDbCommand("SELECT * FROM DNA", m_accessConn);

                // write an sql query into a string. 
                // this query will get all the layers for problem n. 
                //writemessage("query is " + selectCommand.CommandText);
                // execute the query. 
                m_accessConn.Open();


                // Find number of problems
                problemCount = (int)countSequencesCommand.ExecuteScalar();

                OleDbDataReader reader = selectCommand.ExecuteReader();

                if (reader.Read()) // always do this first
                {
                    name = reader.GetString(1);
                    sequence = reader.GetString(2);
                }


            }
            finally
            {
                m_accessConn.Close();
            }

            return new GeneSequence(name, sequence);
        }


        public void WriteGeneSequence(GeneSequence geneSequence)
        {
            //writemessage("writing to database...");
            try
            {
                // because the sequence can be so long, we need to use parameters
                //string insertCommandString = "INSERT INTO DNA (Name, Sequence) VALUES (?, ?)";
                string insertCommandString = "INSERT INTO DNA (Name, Sequence) VALUES (?, ?)";
                //string insertCommandString = "INSERT INTO DNA VALUES (3, ?, ?)";
                OleDbCommand insertCommand = new OleDbCommand(insertCommandString, m_accessConn);
                ASCIIEncoding encoding = new ASCIIEncoding();
                insertCommand.Parameters.Add(new OleDbParameter("name", geneSequence.Name));
                //insertCommand.Parameters.Add(new OleDbParameter("name", encoding.GetBytes(geneSequence.Name.ToCharArray())));
                int a = geneSequence.Name.Length;
                //insertCommand.Parameters.Add(new OleDbParameter("sequence", geneSequence.Sequence.ToCharArray()));
                insertCommand.Parameters.Add(new OleDbParameter("sequence", geneSequence.Sequence));
                m_accessConn.Open();
                insertCommand.ExecuteNonQuery();
            }
            /*catch (Exception e)
            {
                //writemessage("Error trying to write the results to the database");
                //writemessage(e.ToString());
                return;
            }*/
            finally
            {
                m_accessConn.Close();
            }
            //writemessage("done writing to database.  See the row in the tProblems table with problem = " + currentProblem + " to see what happened.");
        }

        public void ConvertTextToDataBaseRows(string pathOfTextFiles)
        {
            foreach(string filename in Directory.GetFiles(pathOfTextFiles))
            {
                System.Console.Out.WriteLine("Reading file: " + filename);
                string name = "";
                string sequence = "";

                // Create an instance of StreamReader to read from a file.
                // The using statement also closes the StreamReader.
                using(StreamReader sr = new StreamReader(filename))
                {
                    name = sr.ReadLine();

                    string line;
                    while((line = sr.ReadLine()) != null)
                    {
                        sequence += line;
                    }
                }
                System.Console.Out.WriteLine("Sequence Length: " + sequence.Length);
                WriteGeneSequence(new GeneSequence(name, sequence));
                //WriteGeneSequence(new GeneSequence("a", sequence));
                System.Console.Out.WriteLine("Wrote GeneSequence.");
            }
        }

        public GeneSequence[] ReadGeneSequences(int max)
        {
            GeneSequence[] result;

            try
            {
                m_accessConn.Open();

                // Find number of problems
                OleDbCommand countSequencesCommand = new OleDbCommand("SELECT MAX(ID) FROM DNA", m_accessConn);
                int sequenceCount = (int)countSequencesCommand.ExecuteScalar();

                if (sequenceCount < max)
                    result = new GeneSequence[sequenceCount];
                else
                    result = new GeneSequence[max];

                // TODO: LIMIT the number of entries returned
                OleDbCommand selectCommand = new OleDbCommand("SELECT TOP " + result.Length + " * FROM DNA ", m_accessConn);
                OleDbDataReader reader = selectCommand.ExecuteReader();

                for (int i = 0; reader.Read() && i < result.Length; ++i)
                {
                    result[i] = new GeneSequence(reader.GetString(1), reader.GetString(2));
                }


            }
            finally
            {
                m_accessConn.Close();
            }

            return result;
        }
    }
}
