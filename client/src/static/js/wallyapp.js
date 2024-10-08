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
const plusButton = document.getElementById('btn-plus')
plusButton.addEventListener('click', zoomIn)
const minusButton = document.getElementById('btn-minus')
minusButton.addEventListener('click', zoomOut)

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


function zoomIn() {
  var url = imgElement.src
  var path = url.replace(new RegExp('/png$'), '').split('-');
  var zoom = Number(path[path.length - 1].replace(new RegExp('^zoom'), ''))
  if (zoom > 0) {  
    zoom = zoom - 1;
  }
  path = path.slice(0, path.length - 1);
  var newUrl = path.join('-') + "-zoom" + zoom.toString() + "/png"    
  const img = new Image();
  img.onload = function() {
      imgElement.width = this.width
      imgElement.height = this.height
  }
  img.src = `${newUrl}`
  imgElement.src = `${newUrl}`
  linkPng.href = `${newUrl}`
}

function zoomOut() {
  var url = imgElement.src
  var path = url.replace(new RegExp('/png$'), '').split('-');
  var zoom = Number(path[path.length - 1].replace(new RegExp('^zoom'), ''))
  if (zoom < 10) {  
    zoom = zoom + 1;
  }
  path = path.slice(0, path.length - 1);
  var newUrl = path.join('-') + "-zoom" + zoom.toString() + "/png"    
  const img = new Image();
  img.onload = function() {
      imgElement.width = this.width
      imgElement.height = this.height
  }
  img.src = `${newUrl}`
  imgElement.src = `${newUrl}`
  linkPng.href = `${newUrl}`    
}

function showExample() {
    var samples = document.getElementById('samples')
    samples.value = 'NA19900\nNA19720\nHG02266\n'
    var regionStart = document.getElementById('regionStart')
    regionStart.value = 33904080
    var regionEnd = document.getElementById('regionEnd')
    regionEnd.value = 33904200
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
