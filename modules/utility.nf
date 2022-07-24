/// This file is only for "utility" processes that are extremely generic.

params.outdir = "output"
params.publishmode = "copy"

process publish_and_rename {

  publishDir params.outdir, mode: params.publishmode
  executor "local"

  input:
    tuple val(filename), path("__infile__")
  output:
    path(filename)

  shell:
    '''
    ln -s $(readlink -f "__infile__") "!{filename}"
    '''
}

process publish {
  publishDir params.outdir, mode: params.publishmode

  executor "local"

  input:
    path filename

  output:
    path filename, includeInputs: true
    
  script:
    """
    """
}

process publish_with_meta {
  publishDir params.outdir, mode: params.publishmode
  executor "local"

  input:
    tuple val(meta), path(filename)

  output:
    path filename, includeInputs: true
    
  script:
  """
  """
}
