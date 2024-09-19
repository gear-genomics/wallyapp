import { saveAs } from 'file-saver'

const API_URL = process.env.API_URL

$('#mainTab a').on('click', function(e) {
  e.preventDefault()
  $(this).tab('show')
})

$('[data-toggle="tooltip"]').tooltip()

const resultLink = document.getElementById('link-results')

const submitButton = document.getElementById('btn-submit')
submitButton.addEventListener('click', function() {
  resultLink.click()
  run()
})

const exampleButton = document.getElementById('btn-example')
exampleButton.addEventListener('click', showExample)

const chromosome = document.getElementById('chromosome')
const dataset = document.getElementById('dataset')
const regionStart = document.querySelector('#regionStart')
const regionEnd = document.querySelector('#regionEnd')
const linkPng = document.getElementById('link-png')
const imgElement = document.getElementById('img-png')
const resultContainer = document.getElementById('result-container')
const resultInfo = document.getElementById('result-info')
const resultError = document.getElementById('result-error')
let downloadUrl

function run() {
  const ds = dataset.querySelector('option:checked').value    
  const chr = chromosome.querySelector('option:checked').value
  const rStart = Number.parseInt(regionStart.value, 10)
  const rEnd = Number.parseInt(regionEnd.value, 10)


  const formData = new FormData()
  formData.append('dataset', ds)
  formData.append('chr', chr)
  formData.append('regionStart', rStart)
  formData.append('regionEnd', rEnd)
  formData.append('samples', document.getElementById('samples').value)

  hideElement(resultContainer)
  hideElement(resultError)
  showElement(resultInfo)

  axios
    .post(`${API_URL}/upload`, formData)
    .then(res => {
      if (res.status === 200) {
          handleSuccess(res.data.data)
      }
    })
    .catch(err => {
      let errorMessage = err
      if (err.response) {
        errorMessage = err.response.data.errors
          .map(error => error.title)
          .join('; ')
      }
      hideElement(resultInfo)
      showElement(resultError)
      resultError.querySelector('#error-message').textContent = errorMessage
    })
}

function handleSuccess(data) {
  hideElement(resultInfo)
  hideElement(resultError)
  showElement(resultContainer)

  downloadUrl = data.url
  linkPng.href = `${API_URL}/${downloadUrl}/png`
  const img = new Image();
  img.onload = function() {
      imgElement.width = this.width
      imgElement.height = this.height
  }
  img.src = `${API_URL}/${downloadUrl}/png`
  imgElement.src = `${API_URL}/${downloadUrl}/png`
}

function showExample() {
    var samples = document.getElementById('samples')
    samples.value = 'NA12878\nHG02611\nNA21099\n'
    var regionStart = document.getElementById('regionStart')
    regionStart.value = 127521550
    var regionEnd = document.getElementById('regionEnd')
    regionEnd.value = 127521650
    var selectchr = document.getElementById('chromosome')
    for (var i = 0 ; i < selectchr.options.length ; i++) {
	if (selectchr.options[i].value == 'chr8') {
	    selectchr.selectedIndex = i;
	}
    }
    var selectds = document.getElementById('dataset')
    for (var i = 0 ; i < selectds.options.length ; i++) {
	if (selectds.options[i].value == '1000 Genomes ONT Vienna GRCh38/hg38') {
	    selectds.selectedIndex = i;
	}
    }
}    

function showElement(element) {
  element.classList.remove('d-none')
}

function hideElement(element) {
  element.classList.add('d-none')
}
